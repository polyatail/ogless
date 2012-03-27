# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "3/16/2012"
__version__ = 0.0

from optparse import OptionParser, OptionGroup
import sys
import os
import shelve
import numpy
from multiprocessing import Pool
import subprocess
import gzip
import time
import tempfile
from cStringIO import StringIO

try:
    import bz2file as bz2
except ImportError:
    sys.stderr.write("WARNING: if using pbzip2-compressed files, specify -z pbzip2 or install bz2file\n")
    import bz2

def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <database.db>",
                          version="%prog " + str(__version__))

    parser.add_option("-a",
                      dest="add",
                      type="str",
                      metavar="[genomes.tab]",
                      default=False,
                      help="add these genomes to the database")

    parser.add_option("-d",
                      dest="delete",
                      type="str",
                      action="append",
                      metavar="genome_name",
                      default=[],
                      help="remove this genome from the database")

    parser.add_option("-c",
                      dest="num_procs",
                      type="int",
                      metavar="[1]",
                      default=1,
                      help="number of concurrent processes")

    parser.add_option("-z",
                      dest="zip_prog",
                      type="str",
                      metavar="[none]",
                      default=False,
                      help="decompress BLAST files with this program")

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
        parser.error("file specified with --add does not exist")

    if options.sbh and options.bbh:
        parser.error("--sbh and --bbh are mutually exclusive")
    elif not (options.sbh or options.bbh):
        options.sbh = True

    if options.add and options.delete:
        parser.error("-a and -d are mutually exclusive")

    if options.zip_prog in ("none", "false", "off"):
        options.zip_prog = False

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

def index_fastas(name_to_data):
    for name, data in name_to_data.items():
        sys.stderr.write("  %s\n" % name)

        genome_num = db["genomes"].index(name)

        for line in open(data["faa"], "r"):
            if line.startswith(">"):
                cds = hash(line[1:-1])

                db["cds_to_genome"][cds] = genome_num
                db["cds_counts"][genome_num] += 1

def parse_blast6_SBH_MP(name_to_data):
    # split into as many workers as we have
    l_name_to_data = name_to_data.items()

    batch_size = len(l_name_to_data) / (options.num_procs * 10)
    batch_size = 1 if batch_size < 1 else batch_size

    batch = []

    for i in range(0, len(l_name_to_data), batch_size):
        batch.append(l_name_to_data[i:i+batch_size])

    # then feed batches to the process pool
    pool = Pool(processes=options.num_procs)
    result = pool.map_async(parse_blast6_SBH, batch)

    # add all the resulting numpy matrices together
    matrices = result.get(2592000)

    sys.stderr.write("adding result matrices\n")

    for i in matrices:
        db["hit_matrix"] += i

def parse_blast6_SBH(name_to_data):
    hits = numpy.zeros((len(db["genomes"]), len(db["genomes"])))

    for name, data in name_to_data:
        sys.stderr.write("  %s\n" % name)

        query_genome = db["genomes"].index(name)

        for line in zip_reader(data["blast6"]):
            line_split = line.rstrip("\n").split("\t")
            
            cds1 = hash(line_split[0])
            ref_genome = db["cds_to_genome"][cds1]

            hits[ref_genome][query_genome] += 1

    return hits

def parse_blast6_BBH_MP(name_to_data):
    # split into as many workers as we have
    l_name_to_data = name_to_data.items()

    batch_size = len(l_name_to_data) / (options.num_procs * 10)
    batch_size = 1 if batch_size < 1 else batch_size

    batch = []

    for i in range(0, len(l_name_to_data), batch_size):
        batch.append(l_name_to_data[i:i+batch_size])

    # then feed batches to the process pool
    pool = Pool(processes=options.num_procs)
    result = pool.map_async(parse_blast6_BBH, batch)

    # without this, workers will stay alive and eat memory when finished
    pool.close()

    # wait for the results
    result_tmpfiles = result.get(2592000)

    sys.stderr.write("processing results\n")
    sys.stderr.write("  loading numpy arrays\n")
    result_data = []

    for i in result_tmpfiles:    
        result_data.append(numpy.load(zip_reader(i)))
        os.unlink(i)

    sys.stderr.write("  concatenating\n")
    result_data = numpy.concatenate(result_data)

    sys.stderr.write("  sorting\n")
    result_data.sort(order="cds_hash")

    sys.stderr.write("  searching for BBHs\n")
    bbh_count = numpy.zeros((len(db["genomes"]), len(db["genomes"])))

    last_hit = (0, 0, 0)

    for cds_hash, genome1, genome2 in result_data:
        if cds_hash == last_hit[0] and \
           last_hit[1:3] in ((genome1, genome2), (genome2, genome1)):
            bbh_count[genome1][genome2] += 1
            bbh_count[genome2][genome1] += 1

        last_hit = (cds_hash, genome1, genome2)

    db["hit_matrix"] += bbh_count

def parse_blast6_BBH(name_to_data):
    result_ndarray = numpy.empty(0, dtype=[("cds_hash", numpy.int64),
                                           ("genome1", numpy.uint16),
                                           ("genome2", numpy.uint16)])

    for name, data in name_to_data:
        sys.stderr.write("  %s\n" % name)

        query_genome = db["genomes"].index(name)
        
        hits = numpy.empty(1000000, dtype=[("cds_hash", numpy.int64),
                                           ("genome1", numpy.uint16),
                                           ("genome2", numpy.uint16)])

        hit_count = 0
    
        for line in zip_reader(data["blast6"]):
            line_split = line.rstrip("\n").split("\t", 3)
            
            cds1 = hash(line_split[0])
            cds2 = hash(line_split[1])
            
            ref_genome = db["cds_to_genome"][cds1]

            try:
                hits[hit_count] = (hash(cds1 + cds2), ref_genome, query_genome)
            except IndexError:
                hits.resize(len(hits) + 100000)
                hits[hit_count] = (hash(cds1 + cds2), ref_genome, query_genome)          

            hit_count += 1

        result_ndarray = numpy.concatenate((result_ndarray, hits[:hit_count]))

    tmpfile = tempfile.NamedTemporaryFile("w", delete=False)
    numpy.save(zip_writer(tmpfile.name), result_ndarray)
    
    return tmpfile.name

class zip_reader():
    def __init__(self, in_fname):
        if options.zip_prog:
            self.proc = subprocess.Popen([options.zip_prog, "-cd", in_fname],
                                    stdout=subprocess.PIPE,
                                    stderr=open(os.devnull, "w"))
    
            self.fp = StringIO(self.proc.stdout.read())
        elif in_fname.lower().endswith((".gz", ".gzip", ".z")):
            self.fp = gzip.open(in_fname, "rb")
        elif in_fname.lower().endswith((".bz2", ".bzip2")):
            self.fp = bz2.BZ2File(in_fname, "rb")
        else:
            self.fp = open(in_fname, "rb")

    def __iter__(self):
        return iter([x for x in self.fp])

    def read(self, size):
        return self.fp.read(size)
        
    def seek(self, offset, whence=0):
        return self.fp.seek(offset, whence)

    def __del__(self):
        self.fp.close()
        
        if options.zip_prog:
            while self.proc.poll() == None:
                time.sleep(0.1)

class zip_writer():
    def __init__(self, out_fname):
        if options.zip_prog:
            self.proc = subprocess.Popen([options.zip_prog, "-c"],
                                         stdin=subprocess.PIPE,
                                         stdout=open(out_fname, "wb"))

            self.fp = self.proc.stdin
        elif out_fname.lower().endswith((".gz", ".gzip", ".z")):
            self.fp = gzip.open(out_fname, "wb")
        elif out_fname.lower().endswith((".bz2", ".bzip2")):
            self.fp = bz2.BZ2File(out_fname, "wb")
        else:
            self.fp = open(out_fname, "wb")
        
    def write(self, data):
        self.fp.write(data)

    def __del__(self):
        self.fp.close()
        
        if options.zip_prog:
            while self.proc.poll() == None:
                time.sleep(0.1)

def new_db(out_fname):
    db = shelve.open(out_fname, flag="n", protocol=-1, writeback=True)

    # index_fastas related stuff    
    db["genomes"] = []
    db["cds_counts"] = numpy.zeros(0)
    db["cds_to_genome"] = {}

    # BLAST parsing related stuff
    db["algorithm"] = "SBH" if options.sbh else "BBH"
    db["hit_matrix"] = numpy.zeros((0, 0))
    
    db.sync()
    
    return db

def main(arguments=sys.argv[1:]):
    parse_options(arguments)

    # cheat by making the database global
    global db

    # make a new database or load an existing one
    if not os.path.exists(args[0]):
        db = new_db(args[0])
    else:
        db = shelve.open(args[0], flag="w", writeback=True)

        if options.sbh and db["algorithm"] != "SBH":
            raise ValueError("Cannot use --sbh with database not generated with SBH")
        elif options.bbh and db["algorithm"] != "BBH":
            raise ValueError("Cannot use --sbh with database not generated with SBH")

    # add some data to the db
    if options.add:
        sys.stderr.write("reading input file\n")
        name_to_data = read_input_file(options.add)

        # make sure we're not going to mess this database up first
        for name in name_to_data.keys():
            if name in db["genomes"]:
                raise ValueError("%s already exists in database" % name)

        # add genomes to the database
        for name in name_to_data.keys():
            db["genomes"].append(name)

        db["cds_counts"].resize(len(db["genomes"]))

        db.sync()

        # index FASTAs
        sys.stderr.write("indexing FASTAs\n")
        index_fastas(name_to_data)
        db.sync()
        
        # resize the hit matrix to accomodate more genomes
        db["hit_matrix"].resize((len(db["genomes"]), len(db["genomes"])))
        
        # and finally, parse the blast6 result files
        sys.stderr.write("parsing blast6 results\n")

        if options.sbh:
            parse_blast6_SBH_MP(name_to_data)
        elif options.bbh:
            parse_blast6_BBH_MP(name_to_data)

        sys.stderr.write("syncing database\n")
        db.sync()

    if options.delete:
        # make sure everything is in the database first
        checked_genomes = {}
        
        for genome in options.delete:
            try:
                genome_num = db["genomes"].index(genome)
            except ValueError:
                raise ValueError("%s is not in the database" % genome)

            checked_genomes[genome] = genome_num

        # then get down to business
        sys.stderr.write("deleting genomes\n")

        for genome, genome_num in checked_genomes.items():
            sys.stderr.write("  %s\n" % genome)

            db["genomes"].remove(genome)

            # remove its columns and rows from the hit matrix
            db["hit_matrix"] = numpy.delete(db["hit_matrix"], genome_num, 0)
            db["hit_matrix"] = numpy.delete(db["hit_matrix"], genome_num, 1)

            # remove cds counts
            db["cds_counts"] = numpy.delete(db["cds_counts"], genome_num)

            sys.stderr.write("    removed cols/rows\n")

            # nuke every CDS belonging to this genome
            nuke_count = 0
            
            for cds_hash in db["cds_to_genome"].keys():
                if db["cds_to_genome"][cds_hash] == genome_num:
                    del db["cds_to_genome"][cds_hash]
                    nuke_count += 1

            sys.stderr.write("    removed %s CDS\n" % nuke_count)

if __name__ == "__main__":
    main()
