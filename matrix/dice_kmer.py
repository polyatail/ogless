# -*- coding: utf-8 -*-

import sys
import os
import numpy
import itertools
from multiprocessing import Process, Queue
from threading import Thread
from Queue import Empty

from Bio import SeqIO

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

def genome2kmers_MP_worker(work_queue, result_queue):
    threads = []

    while True:
        data = work_queue.get(block=True)
        
        if data == False:
            break
        else:
            name, in_fname, length = data

        sys.stderr.write("  %s\n" % name)
    
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
    
def genome2kmers_MP(name_to_data, length):
    work_queue = Queue()
    result_queue = Queue()
    
    for name, data in name_to_data.items()[:100]:
        work_queue.put((name, data["faa"], length))
        
    for i in range(30):
        work_queue.put(False)
        
    processes = [Process(target=genome2kmers_MP_worker, args=(work_queue, result_queue)) for i in range(30)]

    for p in processes:
        p.start()
        
    result_data = []
        
    while True:
        try:
            result_data.append(result_queue.get(block=False))
        except Empty:
            pass
        
        if sum([1 for x in processes if x.is_alive()]) == 0:
            break

    for p in processes:
        p.join()

    return result_data

def kmers2jaccard_MP_worker(work_queue, result_queue):
    threads = []

    while True:
        data = work_queue.get(block=True)
        
        if data == False:
            break
        else:
            result1, result2 = data

            set_intersection = numpy.unique(numpy.intersect1d(result1[1], result2[1]))
            set_union = numpy.unique(numpy.concatenate((result1[1], result2[1])))
    
            jaccard = float(len(set_intersection)) / float(len(set_union))
            dice = float(len(set_intersection)) / float(len(set_intersection) + 0.5 * (len(result1[1]) - len(set_intersection) + len(result2[1]) - len(set_intersection)))

            threads.append(Thread(target=result_queue.put, args=((result1[0], result2[0], dice),)))
            threads[-1].start()
        
    for t in threads:
        t.join()

def kmers2jaccard_MP(kmers):
    work_queue = Queue()
    result_queue = Queue()

    genomes = [x[0] for x in kmers]
    
    for x, y in itertools.combinations(kmers, 2):
        work_queue.put((x, y))
        
    for i in range(30):
        work_queue.put(False)
        
    processes = [Process(target=kmers2jaccard_MP_worker, args=(work_queue, result_queue)) for i in range(30)]

    for p in processes:
        p.start()
        
    result_matrix = numpy.zeros((len(kmers), len(kmers)), dtype=numpy.float64)
        
    while True:
        try:
            data = result_queue.get(block=False)

            result_matrix[genomes.index(data[0])][genomes.index(data[1])] = data[2]
        except Empty:
            pass
        
        if sum([1 for x in processes if x.is_alive()]) == 0:
            break

    for p in processes:
        p.join()

    return result_matrix

def main():
    sys.stderr.write("reading input file\n")
    name_to_data = read_input_file(sys.argv[1])
    
    kmers = genome2kmers_MP(name_to_data, 5)  
    jaccards = kmers2jaccard_MP(kmers)
    
    import pdb; pdb.set_trace()

if __name__ == "__main__":
    main()