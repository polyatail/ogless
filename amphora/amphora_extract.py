# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "3/12/2012"
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

HMMER_PATH = "/usr/bin/hmmscan"
HMMER_DB_PATH = os.path.join(os.path.dirname(__file__), "data", "markers.hmm")

def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <seqs.faa> <label>",
                          version="%prog " + str(__version__))

    parser.add_option("-i",
                      dest="input_file",
                      type="str",
                      metavar="[list.txt]",
                      default=False,
                      help="read list of peptide sequence files")

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

    options, args = parser.parse_args()

    if options.input_file:
        if len(args) <> 0:
            parser.error("incorrect number of arguments")

        if not os.path.exists(options.input_file):
            parser.error("file specified with -i does not exist")
    else:
        if len(args) <> 2:
            parser.error("incorrect number of arguments")

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

def hmmer(in_fname):
    """
    Runs HMMER using the database filename HMMER_DB_PATH and query multifasta filename
    'in_fname' and returns a dict mapping query CDS to their best hits in the
    database. Resolves the situation where a CDS maps to multiple marker genes,
    returning only the best hit (lowest e-value) which spans 70% of the query
    sequence.
    
    Input: query amino acid multifasta filename
    Output: {CDS ID (in in_fname): tuple(marker gene name, evalue)}
    """
    
    # run HMMER on in_fname with reasonable parameters
    out_tmpfile = tempfile.NamedTemporaryFile("r")
    
    hmm_proc = subprocess.call([HMMER_PATH,
                               "-Z", "5000",
                               "-E", "1e-3",
                               "--domtblout", out_tmpfile.name,
                               HMMER_DB_PATH,
                               in_fname],
                               stdout=open(os.devnull, "w"),
                               stderr=subprocess.STDOUT)

    # parse the results when the process is finished
    cds_to_hits = defaultdict(list)

    for line in out_tmpfile:
        if line.startswith("#"):
            continue
        
        line_split = line.rstrip("\n").split()
        
        # only consider hits that span 70% of the query CDS
        if (float(line_split[16]) - float(line_split[15])) / \
           float(line_split[5]) < 0.7:
            continue
        
        cds_to_hits[line_split[3]].append((line_split[0], float(line_split[6])))

    # of the possible marker genes a CDS could be, take the lowest e-value hit
    for cds in cds_to_hits.keys():
        sorted_hits = sorted(cds_to_hits[cds],
                             key=lambda x: x[1],
                             reverse=True)

        cds_to_hits[cds] = sorted_hits[0]

    return cds_to_hits

def extract_genes(cds_to_hits, label, in_fname):
    """
    Extracts the best hits for each marker gene in in_fname and returns a dict
    of amino acid SeqRecords. Resolves the situation where multiple CDS map
    to the same marker gene and returns only the best CDS for each such gene
    
    Input: {CDS ID (in in_fname): tuple(marker gene name, evalue)}
    Output: {marker gene name: SeqRecord(best hit CDS)}
    """
    
    # first, reorganize dict so it maps marker genes to CDS names
    gene_to_hits = defaultdict(list)
    
    for cds, hits in cds_to_hits.items():
        gene_to_hits[hits[0]].append((cds, hits[1]))

    # of the possible CDS a marker gene could map to, take the lowest e-value hit
    for gene in gene_to_hits.keys():
        sorted_hits = sorted(gene_to_hits[gene],
                             key=lambda x: x[1],
                             reverse=True)

        gene_to_hits[gene] = sorted_hits[0]

    # then fetch the FASTA for each CDS best matching each marker gene
    fasta_index = SeqIO.index(in_fname, "fasta")

    gene_to_seqrec = {}
    
    for gene, hits in gene_to_hits.items():
        seqrec = fasta_index[hits[0]]
        seqrec.id = gene + "_" + label
        seqrec.name = gene + "_" + label
        seqrec.description = ""
            
        gene_to_seqrec[gene] = seqrec

    return gene_to_seqrec

def hmm_and_extract(x):
    """
    A simple wrapper around hmmer() and extract_genes() that is compatible with
    pool.map_async()
    
    Input: tuple(genome name, amino acid FASTA filename of CDS)
    Output: tuple(genome_name, {marker gene: SeqRecord(best hit CDS)})
    """
    
    sys.stderr.write("running HMMER on %s\n" % x[0])
    gene_to_hits = hmmer(x[1])
    
    sys.stderr.write("extracting hits from %s\n" % x[0])
    seq_recs = extract_genes(gene_to_hits, x[0], x[1])

    return (x[0], seq_recs)
    
def main(arguments=sys.argv[1:]):
    parse_options(arguments)

    # what files are we evaluating?
    if options.input_file:
        name_to_faa = read_input_file(options.input_file)
    else:
        name_to_faa = {args[1]: args[0]}

    # find marker genes in each
    if options.num_procs > 1:
        # let's have children!
        pool = Pool(processes=options.num_procs)
    
        result = pool.map_async(hmm_and_extract, name_to_faa.items())
        name_to_seq_recs = dict(result.get(2592000))
    else:
        # i'm more old-fashioned, no pool for me
        name_to_seq_recs = {}
        
        for x in name_to_faa.items():
            name_to_seq_recs[x[0]] = hmm_and_extract(x)[1]

    # write the extracted marker gene multifastas
    for name, seq_recs in name_to_seq_recs.items():
        with open(os.path.join(options.output_dir, name + ".faa"), "w") as out_fp:
            SeqIO.write(seq_recs.values(), out_fp, "fasta")

    # make a list of the genomes and their CDS multifastas
    if options.input_file:
        with open(os.path.join(options.output_dir, "genomes.list"), "w") as out_fp:
            for name in name_to_faa.keys():
                row = [name, os.path.join(options.output_dir, name + ".faa")]
                
                out_fp.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()