# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "3/27/2012"
__version__ = 0.0

from optparse import OptionParser, OptionGroup
import sys
import os
import numpy
import scipy
import itertools
import shelve
from multiprocessing import Process, Queue, Array
from threading import Thread
from Queue import Empty
from Bio import SeqIO

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

    parser.add_option("--shmem",
                      dest="shmem",
                      action="store_true",
                      default=False,
                      help="store k-mers in shared memory")

    parser.add_option("-l",
                      dest="kmer_length",
                      type="int",
                      metavar="[5]",
                      default=5,
                      help="length of non-overlapping k-mers")

    options, args = parser.parse_args()

    if len(args) <> 1:
        parser.error("incorrect number of arguments")
        
    if options.add and not os.path.exists(options.add):
        parser.error("file specified with --add does not exist")

    if options.add and options.delete:
        parser.error("-a and -d are mutually exclusive")

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

def genome2kmers_worker(work_queue, result_queue):
    threads = []

    while True:
        data = work_queue.get(block=True)
        
        if data == False:
            break
        else:
            name, in_fname, length = data

        sys.stderr.write("\r  %s/%s %-50s" % (db["genomes"].index(name) + 1, len(db["genomes"]), name[:50]))
    
        kmers = numpy.zeros(1000000, dtype=numpy.int64)
        kmer_count = 0
        
        for seq_rec in SeqIO.parse(in_fname, "fasta"):
            seq = str(seq_rec.seq)
            
            for i in range(0, len(seq) - length, length):
                kmers[kmer_count] = hash(seq[i:i + length])
                kmer_count += 1
                    
        threads.append(Thread(target=result_queue.put, args=((name, numpy.unique(kmers)),)))
        threads[-1].start()
        
    for t in threads:
        t.join()
    
def genome2kmers(name_to_data, length):
    global name_to_kmers

    work_queue = Queue()
    result_queue = Queue()
    
    for name, data in name_to_data.items():
        work_queue.put((name, data["faa"], length))
        
    for i in range(options.num_procs):
        work_queue.put(False)
        
    processes = [Process(target=genome2kmers_worker,
                         args=(work_queue, result_queue)) for i in range(options.num_procs)]

    sys.stderr.write("finding k-mers of length=%s in gene calls\n" % (length,))

    for p in processes:
        p.start()
        
    result_data = []
        
    while True:
        try:
            data = result_queue.get(block=False)
            
            result_data.append((data[0], data[1]))
            db["cds_counts"][db["genomes"].index(data[0])] = len(data[1])
        except Empty:
            pass
        
        if sum([1 for x in processes if x.is_alive()]) == 0:
            break

    for p in processes:
        p.join()

    sys.stderr.write("\n")

    name_to_kmers = dict(result_data)

def kmers2int_worker(work_queue, result_queue):
    threads = []

    while True:
        data = work_queue.get(block=True)
        
        if data == False:
            break
        else:
            for name1, name2 in data:
                kmers1 = name_to_kmers[name1]
                kmers2 = name_to_kmers[name2]

                set_int = numpy.intersect1d(kmers1, kmers2)
                
                result_queue.put((name1, name2, len(set_int)))

def kmers2int():
    work_queue = Queue()
    result_queue = Queue()

    genomes = db["genomes"]
    work_to_do = itertools.combinations(genomes, 2)
    num_work_to_do = int(scipy.comb(len(genomes), 2))
    
    while True:
        batch = list(itertools.islice(work_to_do, 1000))
        
        if not batch:
            break

        work_queue.put(batch)
        
    for i in range(options.num_procs):
        work_queue.put(False)
        
    processes = [Process(target=kmers2int_worker,
                         args=(work_queue, result_queue)) for i in range(options.num_procs)]

    for p in processes:
        p.start()

    sys.stderr.write("filling in matrix\n")
    elements_added = 0

    while True:
        try:
            data = result_queue.get(block=False)
            
            db["hit_matrix"][genomes.index(data[0])][genomes.index(data[1])] = data[2]

            elements_added += 1
        except Empty:
            pass
        
        if sum([1 for x in processes if x.is_alive()]) == 0:
            break

        sys.stderr.write("\r  %s/%s completed" % (elements_added, num_work_to_do))

    for p in processes:
        p.join()
        
    sys.stderr.write("\n")

def mv2shmem():
    sys.stderr.write("copying to shared memory\n")
    
    num_work_to_do = len(name_to_kmers.keys())
    
    for num, name in enumerate(name_to_kmers.keys()):
        name_to_kmers[name] = numpy.frombuffer(Array("l", name_to_kmers[name]).get_obj())

        sys.stderr.write("\r  %s/%s" % (num, num_work_to_do))

    sys.stderr.write("\r  %s/%s\n" % (num, num_work_to_do))

def new_db(out_fname):
    db = shelve.open(out_fname, flag="n", protocol=-1, writeback=True)

    # index_fastas related stuff    
    db["genomes"] = []
    db["cds_counts"] = numpy.zeros(0)

    # hit parsing related stuff
    db["algorithm"] = "kmers"
    db["hit_matrix"] = numpy.zeros((0, 0))
    
    db.sync()
    
    return db

def main(arguments=sys.argv[1:]):
    parse_options(arguments)

    # cheat by making the database and shared memory dict global
    global db, name_to_kmers

    # make a new database or load an existing one
    if not os.path.exists(args[0]):
        db = new_db(args[0])
    else:
        db = shelve.open(args[0], flag="w", writeback=True)
        
        if db["algorithm"] != "kmers":
            raise ValueError("Cannot add to database: not generated with k-mer algorithm")

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

        # resize the hit matrix to accomodate more genomes
        db["hit_matrix"].resize((len(db["genomes"]), len(db["genomes"])))

        genome2kmers(name_to_data, options.kmer_length)
        
        if options.shmem:
            mv2shmem()

        kmers2int()

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


if __name__ == "__main__":
    main()
