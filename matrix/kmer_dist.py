# -*- coding: utf-8 -*-

import sys
import os
import numpy
import scipy
import itertools
from multiprocessing import Process, Queue, Array
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

def genome2kmers_worker(work_queue, result_queue):
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
    
def genome2kmers(name_to_data, length):
    work_queue = Queue()
    result_queue = Queue()
    
    for name, data in name_to_data.items()[:500]:
        work_queue.put((name, data["faa"], length))
        
    for i in range(30):
        work_queue.put(False)
        
    processes = [Process(target=genome2kmers_worker,
                         args=(work_queue, result_queue)) for i in range(30)]

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

    return dict(result_data)

def kmers2distance_worker(work_queue, result_queue):
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
                set_union = numpy.union1d(kmers1, kmers2)

                jaccard = float(len(set_int)) / float(len(set_union))
                dice = float(len(set_int) * 2) / float(len(kmers1) + len(kmers2) - len(set_int) * 2)
    
                result_queue.put((name1, name2, dice, jaccard))

def kmers2distance():
    work_queue = Queue()
    result_queue = Queue()

    genomes = name_to_kmers.keys()
    work_to_do = itertools.combinations(genomes, 2)
    num_work_to_do = scipy.comb(len(genomes), 2)
    
    while True:
        batch = list(itertools.islice(work_to_do, 1000))
        
        if not batch:
            break

        work_queue.put(batch)
        
    for i in range(30):
        work_queue.put(False)
        
    processes = [Process(target=kmers2distance_worker,
                         args=(work_queue, result_queue)) for i in range(30)]

    for p in processes:
        p.start()
        
    dice_matrix = numpy.zeros((len(genomes), len(genomes)), dtype=numpy.float64)
    jaccard_matrix = numpy.zeros((len(genomes), len(genomes)), dtype=numpy.float64)

    elements_added = 0

    while True:
        try:
            data = result_queue.get(block=False)

            dice_matrix[genomes.index(data[0])][genomes.index(data[1])] = data[2]
            jaccard_matrix[genomes.index(data[0])][genomes.index(data[1])] = data[3]

            elements_added += 1
        except Empty:
            pass
        
        if sum([1 for x in processes if x.is_alive()]) == 0:
            break

        sys.stderr.write("\r%s/%s completed" % (elements_added, num_work_to_do))

    for p in processes:
        p.join()

    return dice_matrix, jaccard_matrix

def mv2shmem():
    for name in name_to_kmers.keys():
        name_to_kmers[name] = numpy.frombuffer(Array("l", name_to_kmers[name]).get_obj())

def main():
    sys.stderr.write("reading input file\n")
    name_to_data = read_input_file(sys.argv[1])
    
    global name_to_kmers    
    
    name_to_kmers = genome2kmers(name_to_data, 5)
    
    mv2shmem()

    jaccards = kmers2distance()
    
    import pdb; pdb.set_trace()

if __name__ == "__main__":
    main()