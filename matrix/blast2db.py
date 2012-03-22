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
import gzip

try:
    import bz2file as bz2
except ImportError:
    sys.stderr.write("WARNING: if using pbzip2-compressed files, specify -z pbzip2 or install bz2file\n")
    import bz2

def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <database.db>",
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

    parser.add_option("--sbh",
                      dest="sbh",
                      action="store_true",
                      default=False,
                      help="use single-best hit (SBH) algorithm (default)")

    parser.add_option("--bbh",
                      dest="bbh",
                      action="store_true",
                      default=False,
                      help="use bi-directional best-hit (BBH) algorithm")

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

    if options.sbh and options.bbh:
        parser.error("--sbh and --bbh are mutually exclusive")
    elif not (options.sbh or options.bbh):
        options.sbh = True

def read_input_file(in_fname):
    name_to_data = {}

    for line_num, line in enumerate(open(in_fname, "r")):
        if line.startswith("#"):
            continue

        line_split = line.rstrip("\n").split("\t")
        
        if len(line_split) <> 3:
            raise ValueError("%s:%s: invalid number of fields" % (in_fname, line_num + 1))

        if not os.path.exists(line_split[1]):
            raise ValueError("%s:%s: file does not exist" % (in_fname, line_num + 1))

        if not os.path.exists(line_split[2]):
            raise ValueError("%s:%s: file does not exist" % (in_fname, line_num + 1))

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

def parse_blast6_SBH_MP(name_to_data, db):
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
    result = pool.map_async(parse_blast6_SBH, batch)

    # add all the resulting numpy matrices together
    matrices = result.get(2592000)

    sys.stderr.write("adding result matrices\n")

    for i in matrices:
        db["blast_hits"] += i

def parse_blast6_SBH(in_tuple):
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

def parse_blast6_BBH_MP(name_to_data, db):
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
    result = pool.map_async(parse_blast6_BBH, batch)

    # without this, workers will stay alive and eat memory when finished
    pool.close()

    # wait for the results
    result_data = result.get(2592000)

    sys.stderr.write("processing results\n")
    sys.stderr.write("  concatenating\n")
    result_data = numpy.concatenate(result_data)

    sys.stderr.write("  sorting\n")
    result_data.sort(order="cds_hash")

    sys.stderr.write("  searching for BBHs\n")
    bbh_count = numpy.zeros((len(db["genomes"]), len(db["genomes"])))

    last_hit = (0, 0, 0)

    for x, y, z in result_data:
        # BBH
        if x == last_hit[0] and \
           last_hit[1] in (y, z) and \
           last_hit[2] in (y, z):
            bbh_count[y][z] += 1
            bbh_count[z][y] += 1

        last_hit = (x, y, z)

    db["blast_hits"] += bbh_count

def parse_blast6_BBH(in_tuple):
    name_to_data, cds_to_genome, genome_to_num  = in_tuple 

    result = numpy.empty(0,
                         dtype=[("cds_hash", numpy.int64),
                                ("genome1", numpy.uint16),
                                ("genome2", numpy.uint16)])

    for name, data in name_to_data:
        sys.stderr.write("  %s\n" % name)

        query_genome = genome_to_num[name]
        
        hits = numpy.empty(1000000,
                           dtype=[("cds_hash", numpy.int64),
                                  ("genome1", numpy.uint16),
                                  ("genome2", numpy.uint16)])

        hit_count = 0
    
        for line in zip_wrapper(data["blast6"]):
            line_split = line.rstrip("\n").split("\t", 3)
            
            cds1 = hash(line_split[0])
            cds2 = hash(line_split[1])
            
            try:
                ref_genome = cds_to_genome[cds1]
            except KeyError:
                continue

            try:
                hits[hit_count] = (hash(cds1 + cds2), ref_genome, query_genome)
            except IndexError:
                hits.resize(len(hits) + 100000)
                hits[hit_count] = (hash(cds1 + cds2), ref_genome, query_genome)          

            hit_count += 1

        result = numpy.concatenate((result, hits[:hit_count]))

    return result

def zip_wrapper(in_fname):
    if options.zip_prog:
        proc = subprocess.Popen([options.zip_prog, "-cd", in_fname],
                                stdout=subprocess.PIPE,
                                stderr=open(os.devnull, "w"))

        return [x for x in proc.stdout]
    elif in_fname.lower().endswith((".gz", ".gzip", ".z")):
        return [x for x in gzip.open(in_fname, "r")]
    elif in_fname.lower().endswith((".bz2", ".bzip2")):
        return [x for x in bz2.BZ2File(in_fname, "r")]
    else:
        return [x for x in open(in_fname, "r")]

def new_db(out_fname):
    db = shelve.open(out_fname, flag="n", protocol=-1, writeback=True)
    
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
        db = shelve.open(args[0], flag="w", writeback=True)

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

        if options.sbh:
            if options.num_procs > 1:
                parse_blast6_SBH_MP(name_to_data, db)
            else:
                hits = parse_blast6_SBH((name_to_data.items(), db["cds_to_genome"], db["genome_to_num"]))
                db["blast_hits"] += hits
        elif options.bbh:
            if options.num_procs > 1:
                parse_blast6_BBH_MP(name_to_data, db)
            else:
                hits = parse_blast6_BBH((name_to_data.items(), db["cds_to_genome"], db["genome_to_num"]))
                db["blast_hits"] += hits

        db.sync()

if __name__ == "__main__":
    main()
