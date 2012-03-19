# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "3/13/2012"
__version__ = 0.0

from optparse import OptionParser, OptionGroup
import sys
import os
import subprocess
from collections import defaultdict
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool

MUSCLE_PATH = "/usr/bin/muscle"

def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <genomes.list>",
                          version="%prog " + str(__version__))

    parser.add_option("-o",
                      dest="output_dir",
                      type="str",
                      metavar="[./amphora_out]",
                      default="./amphora_out",
                      help="write output files to this directory")

    parser.add_option("-c",
                      dest="num_procs",
                      type="int",
                      metavar="[1]",
                      default=1,
                      help="number of concurrent processes")

    parser.add_option("--min-markers",
                      dest="min_markers",
                      type="int",
                      metavar="[29]",
                      default=29,
                      help="minimum number of usable markers to run analysis")

    options, args = parser.parse_args()

    if len(args) <> 1:
        parser.error("incorrect number of arguments")

    if not os.path.exists(args[0]):
        parser.error("input file specified does not exist")

    # make the output directory
    options.output_dir = os.path.realpath(options.output_dir)

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

def read_input_file(in_fname):
    name_to_faa = {}

    for line_num, line in enumerate(open(in_fname, "r")):
        if line.startswith("#"):
            continue

        line_split = line.rstrip("\n").split("\t")
        
        if len(line_split) <> 2:
            raise ValueError("%s:%s: invalid number of fields" % (in_fname, line_num+1))

        if not os.path.exists(line_split[1]):
            raise ValueError("%s:%s: file does not exist" % (in_fname, line_num+1))

        name_to_faa[line_split[0]] = line_split[1]

    return name_to_faa

def suggest_exclusions(name_to_seq_recs):
    name_to_genes = {}
    
    for name, seq_recs in name_to_seq_recs.items():
        name_to_genes[name] = set(seq_recs.keys())

    markers = set.intersection(*name_to_genes.values())
    
    if len(markers) >= options.min_markers:
        sys.stderr.write("not necessary to exclude genomes, count %s > %s" % (len(markers), options.min_makers))
        return name_to_genes.keys()

    name_to_score = defaultdict(int)
    
    for nameA in name_to_genes.keys():
        for nameB in set(name_to_genes.keys()).difference([nameA]):
            genes_in_common = name_to_genes[nameA].intersection(name_to_genes[nameB])

            name_to_score[nameA] += len(genes_in_common)
            name_to_score[nameB] += len(genes_in_common)

    # only print the best improvements
    name_to_score = sorted(name_to_score.items(), key=lambda x: x[1])

    tmp_name_to_genes = name_to_genes.copy()

    for name in [x for x, y in name_to_score]:
        del tmp_name_to_genes[name]
        markers = set.intersection(*tmp_name_to_genes.values())
        
        if len(markers) >= options.min_markers:
            sys.stderr.write("removing these genomes increases marker count to %s:\n" % (len(markers),))

            for genome in set(name_to_genes.keys()).difference(tmp_name_to_genes.keys()):
                sys.stderr.write("  " + genome + "\n")

            break

    return tmp_name_to_genes.keys(), markers

def muscle(in_fname):
    # run MUSCLE on in_fname with default parameters
    out_tmpfile = tempfile.NamedTemporaryFile("r")
    
    hmm_proc = subprocess.call([MUSCLE_PATH,
                                "-in", in_fname,
                                "-out", out_tmpfile.name],
                               stdout=open(os.devnull, "w"),
                               stderr=subprocess.STDOUT)

    # return the result as a list of SeqRecords
    return [x for x in SeqIO.parse(out_tmpfile.name, "fasta")]
    
def align_marker(marker_name, name_to_seq_recs):
    sys.stderr.write("  %s\n" % marker_name)

    # pull this marker out of every genome
    query_seq_recs = []
    
    for name, index in name_to_seq_recs.items():
        seq_rec = index[marker_name]
        
        # make sure it's labeled so we can find it later
        seq_rec.id = marker_name + "_" + name
        seq_rec.name = marker_name + "_" + name
        seq_rec.description = ""
        
        query_seq_recs.append(seq_rec)

    # dump sequences to a temporary file
    muscle_in_tempfile = tempfile.NamedTemporaryFile("w")
    SeqIO.write(query_seq_recs, muscle_in_tempfile, "fasta")
    muscle_in_tempfile.flush()
    
    # run muscle
    muscle_out_seq_recs = muscle(muscle_in_tempfile.name)

    name_to_marker = defaultdict(str)
    
    for seq_rec in muscle_out_seq_recs:
        seq_rec_name = "_".join(seq_rec.id.split("_")[1:])
        
        if seq_rec_name not in name_to_seq_recs:
            raise ValueError("%s in MUSCLE output but not in input" % seq_rec_name)

        name_to_marker[seq_rec_name] = str(seq_rec.seq)

    return name_to_marker

def align_marker_MP(x):
    return (x[0], align_marker(*x))

def main(arguments=sys.argv[1:]):
    parse_options(arguments)

    # what files are we evaluating?
    name_to_faa = read_input_file(options.input_file)

    # load name -> seq_recs dict
    name_to_seq_recs = {}
    
    for name, faa in name_to_faa.items():
        name_to_seq_recs[name] = {}

        for seq_rec in SeqIO.parse(faa, "fasta"):
            seq_rec_name = seq_rec.id.split("_")[0]
            
            name_to_seq_recs[name][seq_rec_name] = seq_rec

    # we might need to exclude some genomes with too few markers, and some
    # markers won't appear in all genomes and likewise must be excluded
    kept_names, kept_markers = suggest_exclusions(name_to_seq_recs)

    for name in name_to_seq_recs.keys():
        if name not in kept_names:
            del name_to_seq_recs[name]

    # align each marker individually
    sys.stderr.write("aligning marker genes\n")

    if options.num_procs > 1:
        # let's have children!
        pool = Pool(processes=options.num_procs)
    
        result = pool.map_async(align_marker_MP,
                                [(x, name_to_seq_recs) for x in kept_markers])
        marker_to_seqs = dict(result.get(2592000))
    else:
        marker_to_seqs = {}
        
        for marker in kept_markers:
            marker_to_seqs[marker] = align_marker(marker, name_to_seq_recs)
    
    # combine all marker genes together for each genome
    name_to_cat_seq_recs = defaultdict(str)
    
    for marker in kept_markers:
        for name in kept_names:
            name_to_cat_seq_recs[name] += marker_to_seqs[marker][name]

    for name in kept_names:
        seq_rec = SeqRecord(id=name,
                            description="",
                            seq=Seq(name_to_cat_seq_recs[name]))
                            
        name_to_cat_seq_recs[name] = seq_rec

    with open(os.path.join(options.output_dir, "concat_markers.faa"), "w") as out_fp:
        SeqIO.write(name_to_cat_seq_recs.values(), out_fp, "fasta")

if __name__ == "__main__":
    main()    
    