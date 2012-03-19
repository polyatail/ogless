# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "3/16/2012"
__version__ = 0.0

from optparse import OptionParser, OptionGroup
import sys
import os
from collections import defaultdict
import shelve
import itertools
import numpy
from multiprocessing import Pool
import subprocess
import bz2file
import gzip

def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <tversky.db>",
                          version="%prog " + str(__version__))

    parser.add_option("--add",
                      dest="add",
                      type="str",
                      metavar="[genomes.tab]",
                      default=False,
                      help="add these genomes to the database")

    parser.add_option("--new",
                      dest="new",
                      action="store_true",
                      default=False,
                      help="create a new database")

    parser.add_option("-f",
                      dest="force",
                      action="store_true",
                      default=False,
                      help="force overwrite")

    parser.add_option("-c",
                      dest="num_procs",
                      type="int",
                      metavar="[1]",
                      default=1,
                      help="number of concurrent processes")

    parser.add_option("-z",
                      dest="zip_prog",
                      type="str",
                      metavar="[pbzip2]",
                      default=False,
                      help="decompress blast files with this program")

    options, args = parser.parse_args()

    if len(args) <> 1:
        parser.error("incorrect number of arguments")
        
    if options.add and not os.path.exists(options.add):
        parser.error("--add file does not exist")

    if options.new:
        if os.path.exists(args[0]) and not options.force:
            parser.error("--new specified but database already exists, use -f to overwrite")
    elif not os.path.exists(args[0]):
        parser.error("database file does not exist")

def read_input_file(in_fname):
    name_to_data = {}

    for line_num, line in enumerate(open(in_fname, "r")):
        if line.startswith("#"):
            continue

        line_split = line.rstrip("\n").split("\t")
        
        if len(line_split) <> 3:
            raise ValueError("%s:%s: invalid number of fields" % (in_fname, line_num+1))

        if not os.path.exists(line_split[1]):
            raise ValueError("%s:%s: file does not exist" % (in_fname, line_num+1))

        if not os.path.exists(line_split[2]):
            raise ValueError("%s:%s: file does not exist" % (in_fname, line_num+1))

        name_to_data[line_split[0]] = {"faa": line_split[1],
                                       "blast6": line_split[2]}

    return name_to_data

def index_fastas(name_to_data, db):
    cds_to_genome = db["cds_to_genome"]
    genome_to_cds = db["genome_to_cds"]
    genome_to_num = db["genome_to_num"]
    
    for name, data in name_to_data.items():
        sys.stderr.write("  %s\n" % name)

        genome_num = genome_to_num[name]

        for line in open(data["faa"], "r"):
            if line.startswith(">"):
                cds = hash(line[1:-1])

                cds_to_genome[cds] = genome_num
                genome_to_cds[genome_num].append(cds)

def parse_blast6_MP(name_to_data, db):
    # split into as many workers as we have
    l_name_to_data = name_to_data.items()
    batch_size = len(l_name_to_data) / options.num_procs

    batch = []

    for i in range(0, len(l_name_to_data), batch_size):
        batch.append((l_name_to_data[i:i+batch_size],
                      db["cds_to_genome"],
                      db["genome_to_num"]))

    # then feed batches to the process pool
    pool = Pool(processes=options.num_procs)
    result = pool.map_async(parse_blast6_MP_worker, batch)

    # add all the resulting numpy matrices together
    for i in result.get(2592000):
        db["blast_hits"] += i

def parse_blast6_MP_worker(in_tuple):
    name_to_data, cds_to_genome, genome_to_num  = in_tuple 

    hits = numpy.zeros((len(genome_to_num), len(genome_to_num)))

    for name, data in name_to_data:
        sys.stderr.write("  %s\n" % name)

        query_genome = genome_to_num[name]

        for line in zip_wrapper(data["blast6"]):
            line_split = line.rstrip("\n").split("\t")
            
            cds1 = hash(line_split[0])
            ref_genome = cds_to_genome[cds1]

            hits[ref_genome][query_genome] += 1

    return hits

def zip_wrapper(in_fname):
    if options.zip_prog:
        proc = subprocess.Popen([options.zip, "-cd", in_fname],
                                stdout=subprocess.PIPE,
                                stderr=open(os.devnull, "w"))

        for line in proc.stdout:
            yield line
    elif in_fname.lower().endswith((".gz", ".gzip", ".z")):
        for line in gzip.open(in_fname, "r"):
            yield line
    elif in_fname.lower().endswith((".bz2", ".bzip2")):
        for line in bz2file.BZ2File(in_fname, "r"):
            yield line
    else:
        for line in open(in_fname, "r"):
            yield line

def new_db(out_fname):
    db = shelve.open(out_fname, flag="n", protocol=1, writeback=True)
    
    db["genomes"] = []
    db["genome_to_num"] = {}

    db["cds_to_genome"] = {}
    db["genome_to_cds"] = defaultdict(list)

    db["blast_hits"] = numpy.zeros((0, 0))
    
    db.sync()
    
    return db

def main(arguments=sys.argv[1:]):
    parse_options(arguments)

    # make a new database or load an existing one
    if options.new:
        db = new_db(args[0])
    else:
        db = shelve.open(args[0], flag="w", protocol=1, writeback=True)

    # add some data to the db
    if options.add:
        sys.stderr.write("reading input file\n")
        name_to_data = read_input_file(options.add)

        # make sure we're not going to mess this database up first
        for name in name_to_data.keys():
            if name in db["genomes"]:
                raise ValueError("%s already exists in database" % name)

        # then start making changes
        for name in name_to_data.keys():
            db["genomes"].append(name)

        db["genome_to_num"] = dict([(y, x) for x, y in enumerate(db["genomes"])])

        db.sync()

        sys.stderr.write("indexing FASTAs\n")
        index_fastas(name_to_data, db)
        db.sync()
        
        db["blast_hits"].resize((len(db["genomes"]), len(db["genomes"])))
        
        sys.stderr.write("parsing blast6 results\n")

        if options.num_procs > 1:
            parse_blast6_MP(name_to_data, db)
        else:
            hits = parse_blast6_MP_worker(name_to_data, db["cds_to_genome"], db["genome_to_num"])
            db["blast_hits"] += hits

        db.sync()

if __name__ == "__main__":
    main()
