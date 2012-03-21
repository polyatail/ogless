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

    parser.add_option("--legend",
                      dest="legend",
                      action="store_true",
                      default=False,
                      help="draw a legend for color highlights")

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

    parser.add_option("--highres",
                      dest="highres",
                      action="store_true",
                      default=False,
                      help="draw a large high-res graph (10000x10000p)")

    options, args = parser.parse_args()

    if len(args) <> 1:
        parser.error("incorrect number of arguments")
        
    if not os.path.exists(args[0]):
        parser.error("matrix file doest not exist")

    if options.filter > 1 or \
       options.filter < 0:
        parser.error("filter must be between 0 and 1")

    if options.legend and not options.highlight:
        parser.error("--legend requires --highlight")

    # make the output directory
    options.output_dir = os.path.realpath(options.output_dir)

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

def read_input_file(in_fname, num_fields):
    name_to_data = {}

    for line_num, line in enumerate(open(in_fname, "r")):
        if line.startswith("#") or \
           not line.strip():
            continue

        line_split = line.rstrip("\n").split("\t")
        
        if len(line_split) <> num_fields:
            raise ValueError("%s:%s: invalid number of fields" % (in_fname, line_num + 1))

        name_to_data[line_split[0]] = line_split[1:num_fields]

    return name_to_data

def graph_from_matrix(iterable):
    graph = pygraphviz.AGraph()

    graph.graph_attr["overlap"] = "false"
    graph.graph_attr["forcelabels"] = "false"
    graph.graph_attr["model"] = "mds"
    graph.graph_attr["outputorder"] = "edgesfirst"
    graph.graph_attr["packmode"] = "graph"

    graph.edge_attr["penwidth"] = "0.5"
    graph.edge_attr["color"] = "#00000030"
    
    graph.node_attr["style"] = "filled"
    graph.node_attr["color"] = "#0000005f"
    graph.node_attr["shape"] = "point"
    graph.node_attr["fixedsize"] = "true"
    graph.node_attr["width"] = "0.1"
    graph.node_attr["height"] = "0.1"

    # read header and generate a dict mapping position in header to genome name
    header = iterable.next()
    
    if header.startswith("name"):
        header = header.rstrip("\n").split("\t")[1:]

    num_to_name = dict(enumerate(header))

    # read matrix line-by-line, keeping only edges with scores above options.filter
    sys.stderr.write("loading matrix into graph\n")

    for line_num, line in enumerate(iterable):
        line_split = line.rstrip("\n").split("\t")

        if len(line_split[1:]) <> len(header):
            raise ValueError("%s:%s: invalid number of fields" % (args[0], line_num + 1))

        kept_fields = filter(lambda (x, y): y >= options.filter, enumerate(map(float, line_split[1:])))
        
        for num, score in kept_fields:
            graph.add_edge(line_split[0], num_to_name[num], len=str(1 - score))

        sys.stderr.write("\r  %s%% complete" % round(100 * (float(line_num) / float(len(header)))))

    sys.stderr.write("\r  100% complete\n")

    return graph

def color_nodes(graph, name_to_highlight):
    for node in graph.nodes():
        for regexp, (name, color) in name_to_highlight.items():
            if re.match(regexp, node.name):
                node.attr["color"] = color
                break
        else:
            node.attr["color"] = name_to_highlight["default"][1]

def label_nodes(graph):
    for node in graph.nodes():
        node.attr["xlabel"] = node.name

def print_matrix(db, matrix):
    for genome_num, genome_name in enumerate(db["genomes"]):
        row = [genome_name] + [x[:8] for x in map(str, matrix[genome_num])]
        print "\t".join(row)

def make_legend(graph, name_to_highlight):
    table = \
"""
<table border='0' cellborder='1' cellspacing='0' cellpadding='4' color='#000000' bgcolor='#FFFFFF'>
  <tr>
    <td colspan='2' valign='middle'><b><font face='Helvetica'>Legend</font></b></td>
  </tr>
  %s
</table>
"""

    entry = \
"""
  <tr>
    <td bgcolor='%s'>&nbsp;</td>
    <td valign='middle'><font face='Helvetica'>%s</font></td>
  </tr>
"""

    items = [entry % (name_to_highlight["default"][1], name_to_highlight["default"][0])]
   
    for regexp, (name, color) in sorted(name_to_highlight.items(), key=lambda x: x[1][0]):
        if regexp == "default":
            continue
        
        items.append(entry % (color, name))

    label = (table % "".join(items)).replace("\n", "")

    graph.node_attr["fixedsize"] = "false"
    graph.add_node("legend", color="#00000000", shape="none", margin="0.5", label="<" + label + ">")

def main(arguments=sys.argv[1:]):
    parse_options(arguments)

    graph = graph_from_matrix(open(args[0], "r"))

    if options.highlight:
        sys.stderr.write("adding highlight colors\n")
        name_to_highlight = read_input_file(options.highlight, 3)

        color_nodes(graph, name_to_highlight)

        if options.legend:
            make_legend(graph, name_to_highlight)

    if options.label:
        sys.stderr.write("adding node labels\n")
        label_nodes(graph)

    if options.highres:
        graph.edge_attr["pen_width"] = "2"
        
        graph.graph_attr["size"] = "'100,100!'"
        graph.graph_attr["ratio"] = "fill"
        graph.graph_attr["sep"] = "'+15,15'"

        graph.node_attr["height"] = "0.4"
        graph.node_attr["width"] = "0.4"

    graph.write(os.path.join(options.output_dir, "graph.dot"))
    
    if not options.nodraw:
        graph.draw(os.path.join(options.output_dir, "graph.png"), format="png", prog="neato")

if __name__ == "__main__":
    main()