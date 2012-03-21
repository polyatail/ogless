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

def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <tversky.db>",
                          version="%prog " + str(__version__))

    parser.add_option("-a",
                      dest="alpha",
                      type="float",
                      metavar="[0.5]",
                      default=0.5,
                      help="alpha value for Tversky calculation")

    options, args = parser.parse_args()

    if len(args) <> 1:
        parser.error("incorrect number of arguments")
        
    if not os.path.exists(args[0]):
        parser.error("database file does not exist")

    if options.alpha > 1 or \
       options.alpha < 0:
        parser.error("alpha must be between 0 and 1")

def tversky_matrix(db, alpha):
    tversky_matrix = numpy.zeros((len(db["genomes"]), len(db["genomes"])))

    for nameA, nameB in itertools.combinations(db["genomes"], 2):
        nameA_num = db["genome_to_num"][nameA]
        nameB_num = db["genome_to_num"][nameB]
        
        sizeA = len(db["genome_to_cds"][nameA_num])
        sizeB = len(db["genome_to_cds"][nameB_num])

        hitsA = db["blast_hits"][nameA_num][nameB_num]
        hitsB = db["blast_hits"][nameB_num][nameA_num]

        if hitsA > sizeA or \
           hitsB > sizeB:
            raise ValueError("shared genes cannot exceed genome size")

        shared_genes = min([hitsA, hitsB])

        a_not_b = sizeA - shared_genes
        b_not_a = sizeB - shared_genes
        
        tverskyAB = float(shared_genes) / float(shared_genes + alpha * a_not_b + (1 - alpha) * b_not_a)
        tverskyBA = float(shared_genes) / float(shared_genes + (1 - alpha) * a_not_b + alpha * b_not_a)

        tversky_matrix[nameA_num][nameB_num] = tverskyAB
        tversky_matrix[nameB_num][nameA_num] = tverskyBA

    return tversky_matrix

def print_matrix(db, matrix):
    print "\t".join(["name"] + db["genomes"])

    for genome_num, genome_name in enumerate(db["genomes"]):
        row = [genome_name] + [x[:8] for x in map(str, matrix[genome_num])]
        print "\t".join(row)

def main(arguments=sys.argv[1:]):
    parse_options(arguments)

    db = shelve.open(args[0], flag="w", protocol=1, writeback=True)

    matrix = tversky_matrix(db, options.alpha)
    
    print_matrix(db, matrix)

if __name__ == "__main__":
    main()


