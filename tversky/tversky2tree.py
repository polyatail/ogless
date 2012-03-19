# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "3/16/2012"
__version__ = 0.0

from optparse import OptionParser, OptionGroup
import sys
import os
import subprocess
import tempfile

QUICKTREE_PATH = "/usr/bin/quicktree"

def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <tversky.matrix>",
                          version="%prog " + str(__version__))

    parser.add_option("-o",
                      dest="output_dir",
                      type="str",
                      metavar="[./tversky_out]",
                      default="./tversky_out",
                      help="write output files to this directory")

    options, args = parser.parse_args()

    if len(args) > 1:
        parser.error("incorrect number of arguments")
    elif len(args) == 1 and not os.path.exists(args[0]):
        parser.error("matrix file does not exist")

def quicktree_matrix(iterable):
    # tempfile stores matrix in PHYLIP format
    tmp_out = tempfile.NamedTemporaryFile("w")

    # read the first line (header)
    header = iterable.next()
    
    if header.startswith("name"):
        header = header.rstrip("\n").split("\t")[1:]

    tmp_out.write("\t%s\n" % len(header))

    # we're going to have to rename every line. quicktree can't handle long
    # species/genome names in the first column
    name_to_num = dict([(y, x) for x, y in enumerate(header)])

    for line_num, line in enumerate(iterable):
        line_split = line.rstrip("\n").split("\t")

        if len(line_split[1:]) <> len(header):
            raise ValueError("%s:%s: invalid number of fields" % (args[0], line_num + 1))

        # rename this column
        line_split[0] = name_to_num[line_split[0]]

        # take the inverse of every element as tversky reflects similarity
        # and we'd like to consider distance
        line_split[1:] = [(1 - float(x)) for x in line_split[1:]]

        tmp_out.write("\t".join(map(str, line_split)) + "\n")

    return tmp_out, name_to_num

def quicktree(in_fname, name_to_num):
    # flip name_to_num so we can undo the line renaming we did above
    num_to_name = dict([(y, x) for x, y in name_to_num.items()])

    proc = subprocess.Popen([QUICKTREE_PATH,
                             "-in", "m",
                             "-out", "t",
                             in_fname],
                            stdout=subprocess.PIPE,
                            stderr=open(os.devnull, "w"))

    # rename the lines as they come off the press
    tree = []

    for line in proc.stdout:
        line_split = line.rstrip("\n").split(":")
        
        if len(line_split) == 2:
            if line_split[0]:
                line_split[0] = num_to_name[int(line_split[0])]

        tree.append(":".join(line_split))

    return "".join(tree)

def main(arguments=sys.argv[1:]):
    parse_options(arguments)

    if len(args) == 1:
        in_fp = open(args[0], "r")
    else:
        in_fp = sys.stdin

    tmp_out, name_to_num = quicktree_matrix(in_fp)

    tree = quicktree(tmp_out.name, name_to_num)

    print tree

if __name__ == "__main__":
    main()