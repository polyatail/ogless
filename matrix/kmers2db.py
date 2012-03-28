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
import time
from multiprocessing import Process, Queue, RawArray
from threading import Thread
from Queue import Empty, Queue as ThreadQueue
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
                      dest="num_tasks",
                      type="int",
                      metavar="[1]",
                      default=1,
                      help="number of concurrent tasks")

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
                      help="length of k-mers")

    parser.add_option("-b",
                      dest="bump",
                      type="int",
                      metavar="[5]",
                      default=5,
                      help="skip this many positions between k-mers")

    options, args = parser.parse_args()

    if len(args) <> 1:
        parser.error("incorrect number of arguments")
        
    if options.add and not os.path.exists(options.add):
        parser.error("file specified with --add does not exist")

    if options.add and options.delete:
        parser.error("-a and -d are mutually exclusive")
        
    if options.kmer_length < 1:
        parser.error("-l must be at least 1")
        
    if options.bump < 1:
        parser.error("-b must be at least 1")

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
            num, in_fname, length = data

        sys.stderr.write("\r%-79s" % ("  %s/%s %-50s" % (num + 1, len(db["genomes"]),
                                                         db["genomes"][num][:50]),))
    
        kmers = numpy.zeros(1000000, dtype=numpy.int64)
        kmer_count = 0
        
        for seq_rec in SeqIO.parse(in_fname, "fasta"):
            seq = str(seq_rec.seq)
            
            for i in range(0, len(seq) - length, options.bump):
                try:
                    kmers[kmer_count] = hash(seq[i:i + length])
                except IndexError:
                    kmers.resize(len(kmers) + 100000)
                    kmers[kmer_count] = hash(seq[i:i + length])

                kmer_count += 1
                    
        threads.append(Thread(target=result_queue.put, args=((num, numpy.unique(kmers)),)))
        threads[-1].start()
        
    for t in threads:
        t.join()
    
def genome2kmers(name_to_data):
    global num_to_kmers

    work_queue = Queue()
    result_queue = Queue()
    
    for num, data in enumerate(name_to_data.values()):
        work_queue.put((num, data["faa"], options.kmer_length))
        
    for i in range(options.num_tasks):
        work_queue.put(False)
        
    tasks = [Process(target=genome2kmers_worker,
                     args=(work_queue, result_queue)) for i in range(options.num_tasks)]

    sys.stderr.write("finding %s-mers skip %s in gene calls\n" % (options.kmer_length,
                                                                  options.bump))

    for t in tasks:
        t.start()
        
    result_data = []
        
    while True:
        try:
            data = result_queue.get(block=False)
            
            db["cds_counts"][data[0]] = len(data[1])
            
            if options.shmem:
                new_array = numpy.frombuffer(RawArray("l", len(data[1])), dtype=numpy.int64)
                new_array[:] = data[1]
                new_array.setflags(write=False)
                
                result_data.append((data[0], new_array))
            else:
                result_data.append((data[0], data[1]))            
        except Empty:
            pass
        
        if sum([1 for x in tasks if x.is_alive()]) == 0:
            break

    for t in tasks:
        t.join()

    sys.stderr.write("\n")

    num_to_kmers = dict(result_data)

def kmers2int_worker(work_queue, result_queue):
    intersect1d = numpy.intersect1d

    while True:
        data = work_queue.get(block=True)
        
        if data == False:
            break
        else:
            for num1, num2 in data:
                int_set = intersect1d(num_to_kmers[num1],
                                      num_to_kmers[num2],
                                      assume_unique=True)

                result_queue.put((num1, num2, int_set.size))

def kmers2int():
    work_queue = Queue()
    result_queue = Queue()

    work_to_do = itertools.combinations(range(len(db["genomes"])), 2)
    num_work_to_do = int(round(scipy.comb(len(db["genomes"]), 2)))
    
    while True:
        batch = list(itertools.islice(work_to_do, 1000))
        
        if not batch:
            break

        work_queue.put(batch)
        
    for i in range(options.num_tasks):
        work_queue.put(False)
        
    tasks = [Process(target=kmers2int_worker,
                     args=(work_queue, result_queue)) for i in range(options.num_tasks)]

    for t in tasks:
        t.start()

    sys.stderr.write("populating matrix\n")

    true_s_time = time.time()
    s_time = time.time()
    status_change = True
    work_done = 0
    delta_time = []

    while True:
        for _ in range(10):
            try:
                data = result_queue.get(block=False)

                db["hit_matrix"][data[0]][data[1]] = data[2]

                work_done += 1
                status_change = True
            except Empty:
                break

        if sum([1 for t in tasks if t.is_alive()]) == 0:
            break

        if status_change and work_done % 100 == 0:
            try:
                work_per_time = float(100 * len(delta_time)) / float(sum(delta_time))
            except ZeroDivisionError:
                work_per_time = 0

            delta_time.append(time.time() - s_time)
            
            try:
                time_to_finish = (num_work_to_do - work_done) / work_per_time
            except ZeroDivisionError:
                time_to_finish = 0
            
            sys.stderr.write("\r%-79s" % ("  %s/%s completed [%.2f/s, %s remaining]" % (work_done,
                             num_work_to_do, work_per_time, formattime(time_to_finish)),))

            if len(delta_time) == 100:
                delta_time = delta_time[1:]                

            s_time = time.time()
            status_change = False

    for t in tasks:
        t.join()
        
    sys.stderr.write("\r%-79s\n" % ("  %s/%s completed in %s" % (work_done,
                     num_work_to_do, formattime(time.time() - true_s_time)),))

def formattime(secs):
    secs = int(secs)

    hours = secs / 3600
    secs -= 3600 * hours
    minutes = secs / 60
    secs -= 60 * minutes

    if hours == 0 and minutes == 0:
        return "%ds" % (secs,)
    elif hours == 0:
        return "%dm %ds" % (minutes, secs)
    else:
        return "%dh %dm %ds" % (hours, minutes, secs)

def new_db(out_fname):
    db = shelve.open(out_fname, flag="n", protocol=-1, writeback=True)

    # index_fastas related stuff    
    db["genomes"] = []
    db["cds_counts"] = numpy.zeros(0, dtype=numpy.int64)

    # hit parsing related stuff
    db["algorithm"] = "kmers"
    db["hit_matrix"] = numpy.zeros((0, 0), dtype=numpy.int64)
    
    db.sync()
    
    return db

def main(arguments=sys.argv[1:]):
    parse_options(arguments)

    # cheat by making the database and shared memory dict global
    global db, num_to_kmers

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

        genome2kmers(name_to_data)

        kmers2int()
        
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

if __name__ == "__main__":
    main()
