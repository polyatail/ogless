# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "3/16/2012"
__version__ = 0.0

from optparse import OptionParser, OptionGroup, IndentedHelpFormatter
import sys
import os
from collections import defaultdict
import shelve
import itertools
import numpy
import textwrap

class BetterFormatter(IndentedHelpFormatter):
    def format_option(self, option):
        result = []

        opts = self.option_strings[option]

        opt_width = self.help_position - self.current_indent - 2 

        if len(opts) > opt_width:
            opts = "%*s%s\n" % (self.current_indent, "", opts)
            indent_first = self.help_position
        else:
            opts = "%*s%-*s  " % (self.current_indent, "", opt_width, opts)
            indent_first = 0 

        result.append(opts)

        if option.help:
            help_text = self.expand_default(option)
            help_text = self.parser.expand_prog_name(help_text)

            help_lines = []
            wrapper = textwrap.TextWrapper(width=self.help_width)

            for p in map(wrapper.wrap, help_text.split('\n')):
                if p:
                    help_lines.extend(p)
                else:
                    help_lines.append("")

            result.append("%*s%s\n" % (indent_first, "", help_lines[0]))
            result.extend(["%*s%s\n" % (self.help_position, "", line) for line in help_lines[1:]])
        elif opts[-1] != "\n":
            result.append("\n")
            
        return "".join(result)

def parse_options(arguments):
    global options, args

    parser = OptionParser(formatter=BetterFormatter(), usage="%prog [options] <database.db>",
                          version="%prog " + str(__version__))

    parser.add_option("-s",
                      dest="sim_func",
                      type="choice",
                      metavar="[jaccard]",
                      choices=["tversky", "jaccard", "dice"],
                      default="jaccard",
                      help="similarity function to compare genomes\n"
                           "    tversky\n"
                           "    jaccard (default)\n"
                           "    dice")

    parser.add_option("-a",
                      dest="alpha",
                      type="float",
                      metavar="[0.5]",
                      default=False,
                      help="alpha value for Tversky index")

    parser.add_option("-b",
                      dest="beta",
                      type="float",
                      metavar="[0.5]",
                      default=False,
                      help="beta value for Tversky index")

    parser.add_option("-r",
                      dest="resolve",
                      type="choice",
                      metavar="[min]",
                      choices=["min", "mean", "max"],
                      default="min",
                      help="when calculating similarity on a single best hit (SBH)\n"
                           "matrix, resolve differences in the number of orthologs\n"
                           "found by taking the:\n"
                           "    min (default)\n"
                           "    mean\n"
                           "    max")

    options, args = parser.parse_args()

    if len(args) <> 1:
        parser.error("incorrect number of arguments")

    if options.sim_func == "tversky":
        if isinstance(options.alpha, bool):
            options.alpha = float(0.5)
        elif options.alpha > 1 or \
             options.alpha < 0:
            parser.error("alpha must be between 0 and 1")

        if isinstance(options.beta, bool):
            options.beta = float(0.5)
        elif options.beta > 1 or \
             options.beta < 0:
            parser.error("beta must be between 0 and 1")
    elif options.alpha or options.beta:
        parser.error("must specify -s tversky to use -a or -b")
        
    if not os.path.exists(args[0]):
        parser.error("database file does not exist")

def tversky(sizeA, sizeB, size_int):
    return float(size_int) / float(size_int + options.alpha * (sizeA - size_int) + options.beta * (sizeB - size_int))

def jaccard(sizeA, sizeB, size_int):
    return float(size_int) / float(sizeA + sizeB - size_int)

def dice(sizeA, sizeB, size_int):
    return float(size_int) / float(size_int + 0.5 * (sizeA - size_int) + 0.5 * (sizeB - size_int))

def similarity_matrix(db, sim_func, check_func):
    matrix = numpy.zeros((len(db["genomes"]), len(db["genomes"])))

    for nameA, nameB in itertools.combinations(db["genomes"], 2):
        nameA_num = db["genomes"].index(nameA)
        nameB_num = db["genomes"].index(nameB)
        
        sizeA = db["cds_counts"][nameA_num]
        sizeB = db["cds_counts"][nameB_num]

        hits_count = check_func((db["hit_matrix"][nameA_num][nameB_num],
                                 db["hit_matrix"][nameB_num][nameA_num]))

        matrix[nameA_num][nameB_num] = sim_func(sizeA, sizeB, hits_count)
        matrix[nameB_num][nameA_num] = sim_func(sizeB, sizeA, hits_count)

    return matrix

def print_matrix(db, matrix):
    print "\t".join(["name"] + db["genomes"])

    for genome_num, genome_name in enumerate(db["genomes"]):
        row = [genome_name] + [x[:8] for x in map(str, matrix[genome_num])]
        print "\t".join(row)

def main(arguments=sys.argv[1:]):
    parse_options(arguments)
    
    db = shelve.open(args[0], flag="w", protocol=1, writeback=True)

    if options.sim_func == "tversky":
        sim_func = tversky
    elif options.sim_func == "jaccard":
        sim_func = jaccard
    elif options.sim_func == "dice":
        sim_func = dice

    if db["algorithm"] == "SBH":
        if options.resolve == "min":
            check_func = min
        elif options.resolve == "mean":
            check_func = lambda x: float(sum(x)) / float(len(x))
        elif options.resolve == "max":
            check_func = max
    elif db["algorithm"] == "BBH":
        def check_func(x):
            if x[0] != x[1]:
                raise ValueError("BBH matrix not symmetrical")
            
            return x[0]
    else:
        raise ValueError("Unknown algorithm type '%s'" % (db["algorithm"],))

    matrix = similarity_matrix(db, sim_func, check_func)
    
    print_matrix(db, matrix)

if __name__ == "__main__":
    main()


