# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "3/15/2012"
__version__ = 0.0

NUCMER_PATH = "/usr/bin/nucmer"
SHOWCOORDS_PATH = "/usr/bin/show-coords"

from optparse import OptionParser, OptionGroup
import sys
import os
import subprocess
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
import itertools

def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <synteny.tab>",
                          version="%prog " + str(__version__))

    parser.add_option("-o",
                      dest="output_dir",
                      type="str",
                      metavar="[./circos_out]",
                      default="./circos_out",
                      help="write output files to this directory")

    options, args = parser.parse_args()

    if len(args) <> 1:
        parser.error("incorrect number of arguments")

    # make the output directory
    options.output_dir = os.path.realpath(options.output_dir)

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

def read_input_file(in_fname):
    name_to_data = {}

    for line_num, line in enumerate(open(in_fname, "r")):
        if line.startswith("#"):
            continue

        line_split = line.rstrip("\n").split("\t")
        
        if len(line_split) <> 3:
            raise ValueError("%s:%s: invalid number of fields" % (in_fname, line_num+1))

        if not os.path.exists(line_split[2]):
            raise ValueError("%s:%s: file does not exist" % (in_fname, line_num+1))

        name_to_data[line_split[0]] = {"faa": line_split[2],
                                       "color": line_split[1]}

    return name_to_data

def fasta2karyotype(name, color, in_fname):
    genome = []
    
    for seq_rec in SeqIO.parse(in_fname, "fasta"):
        genome.append(str(seq_rec.seq))

    genome_size = sum([len(x) for x in genome])

    # write karyotype
    new_file = [" ".join(["chr -", name, name, "0", str(genome_size), color])]

    running_total = 0

    for contig_num, contig in enumerate(genome):
        contig_name = "contig%s" % (contig_num + 1,)

        new_file.append(" ".join(["band", name, contig_name, contig_name,
                                  str(running_total), str(running_total + len(contig)),
                                  color]))

        running_total += len(contig)

    with open(os.path.join(options.output_dir, "karyotype." + name + ".txt"), "w") as out_fp:
        out_fp.write("\n".join(new_file))

    # write FASTA
    seq_rec = SeqRecord(id=name, name=name, description="", seq=Seq("".join(genome)))
    SeqIO.write(seq_rec, os.path.join(options.output_dir, name + ".fasta"), "fasta")

def nucmer(ref_fname, query_fname):
    out_tmpfile = tempfile.NamedTemporaryFile("w+")

    nucmer_proc = subprocess.call([NUCMER_PATH,
                                   "-b", "1000",
                                   "--maxmatch",
                                   "-p", out_tmpfile.name,
                                   ref_fname, query_fname],
                                  stdout=open(os.devnull, "w"),
                                  stderr=subprocess.STDOUT)

    showcoords_proc = subprocess.Popen([SHOWCOORDS_PATH,
                                        "-T",
                                        out_tmpfile.name + ".delta"],
                                       stdout=subprocess.PIPE,
                                       stderr=open(os.devnull, "w"))

    while showcoords_proc.poll() == None:
        time.sleep(1)

    matches = []

    col_order = ["ref_start", "ref_end", "query_start", "query_end",
                 "ref_length", "query_length", "perc_identity",
                 "ref_name", "query_name"]

    for line in showcoords_proc.stdout:
        try:
            int(line[0])
        except ValueError:
            continue
        
        line_split = line.rstrip("\n").split("\t")
        row = dict(zip(col_order, map(int, line_split[:6]) + [float(line_split[6])] + line_split[7:]))
        matches.append(row)

    return matches

def nucmer2links(matches, ref_name, query_name):
    links = []
    
    for link_count, x in enumerate(matches):
        if x["ref_start"] > x["ref_end"]:
            x["ref_start"], x["ref_end"] = x["ref_end"], x["ref_start"]

        if x["query_start"] > x["query_end"]:
            x["query_start"], x["query_end"] = x["query_end"], x["query_start"]
        
        row1 = ["link%s" % (link_count + 1,), x["ref_name"], x["ref_start"], x["ref_end"]]
        row2 = ["link%s" % (link_count + 1,), x["query_name"], x["query_start"], x["query_end"]]

        links.append(" ".join(map(str, row1)))
        links.append(" ".join(map(str, row2)))

    with open(os.path.join(options.output_dir, "links.%s_%s.txt" % (ref_name, query_name)), "w") as out_fp:
        out_fp.write("\n".join(links))

def main(arguments=sys.argv[1:]):
    parse_options(arguments)

    # what files are we evaluating?
    sys.stderr.write("reading input file\n")
    name_to_data = read_input_file(args[0])

    # generate karyotypes and renamed FASTAS
    sys.stderr.write("generating karyotypes\n")
    name_to_seq_recs = {}

    for name, data in name_to_data.items():
        sys.stderr.write("  %s\n" % name)

        seq_recs = fasta2karyotype(name, data["color"], data["faa"])

    # compare and generate links
    sys.stderr.write("running NUCMER\n")
    
    for ref_name, query_name in itertools.combinations(name_to_data.keys(), 2):
        sys.stderr.write("  %s vs %s\n" % (ref_name, query_name))

        ref_file = os.path.join(options.output_dir, ref_name + ".fasta")
        query_file = os.path.join(options.output_dir, query_name + ".fasta")
        
        nucmer_matches = nucmer(ref_file, query_file)
        nucmer2links(nucmer_matches, ref_name, query_name)

if __name__ == "__main__":
    main()    
