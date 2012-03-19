# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "3/16/2012"
__version__ = 0.0

from optparse import OptionParser, OptionGroup
import sys
import os
import pygraphviz
import re

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

    parser.add_option("-f",
                      dest="filter",
                      type="float",
                      metavar="[0.4]",
                      default=0.4,
                      help="filter edges with scores less than this")

    parser.add_option("--highlight",
                      dest="highlight",
                      type="str",
                      metavar="[highlights.txt]",
                      default=False,
                      help="load node color mapping from this file")

    parser.add_option("--label",
                      dest="label",
                      action="store_true",
                      default=False,
                      help="label nodes")

    parser.add_option("--nodraw",
                      dest="nodraw",
                      action="store_true",
                      default=False,
                      help="save graph in .dot format but do not draw")

    options, args = parser.parse_args()

    if len(args) <> 1:
        parser.error("incorrect number of arguments")
        
    if not os.path.exists(args[0]):
        parser.error("matrix file doest not exist")

    if options.filter > 1 or \
       options.filter < 0:
        parser.error("filter must be between 0 and 1")

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
        
        if len(line_split) <> 2:
            raise ValueError("%s:%s: invalid number of fields" % (in_fname, line_num + 1))

        name_to_data[line_split[0]] = line_split[1]

    return name_to_data

def graph_from_matrix(iterable):
    graph = pygraphviz.AGraph()

    graph.graph_attr["overlap"] = "false"
    graph.graph_attr["forcelabels"] = "false"
    graph.graph_attr["model"] = "mds"
    graph.graph_attr["defaultdist"] = "1"

    graph.edge_attr["penwidth"] = "0.5"
    graph.edge_attr["color"] = "#00000010"
    
    graph.node_attr["style"] = "filled"
    graph.node_attr["color"] = "#0000005f"
    graph.node_attr["shape"] = "point"
    graph.node_attr["fixedsize"] = "true"
    graph.node_attr["width"] = "0.1"
    graph.node_attr["height"] = "0.1"

    header = iterable.next()
    
    if header.startswith("name"):
        header = header.rstrip("\n").split("\t")[1:]

    for line_num, line in enumerate(iterable):
        line_split = line.rstrip("\n").split("\t")

        if len(line_split[1:]) <> len(header):
            raise ValueError("%s:%s: invalid number of fields" % (args[0], line_num + 1))

        ref_genome = line_split[0]
        
        for query_genome, score in zip(header, map(float, line_split[1:])):
            if score < options.filter:
                continue
            
            graph.add_edge(ref_genome, query_genome, len=str(1 - score))

    return graph

def color_nodes(graph, name_to_color):
    for node in graph.nodes():
        for name, color in name_to_color.items():
            if re.match(name, node.name):
                node.attr["color"] = color

def label_nodes(graph):
    for node in graph.nodes():
        node.attr["xlabel"] = node.name

def print_matrix(db, matrix):
    for genome_num, genome_name in enumerate(db["genomes"]):
        row = [genome_name] + [x[:8] for x in map(str, matrix[genome_num])]
        print "\t".join(row)

def main(arguments=sys.argv[1:]):
    parse_options(arguments)

    graph = graph_from_matrix(open(args[0], "r"))

    if options.highlight:
        name_to_color = read_input_file(options.highlight)

        color_nodes(graph, name_to_color)

    if options.label:
        label_nodes(graph)

    graph.write(os.path.join(options.output_dir, "graph.dot"))
    
    if not options.nodraw:
        graph.draw(os.path.join(options.output_dir, "graph.png"), format="png", prog="neato")

if __name__ == "__main__":
    main()