#!/usr/bin/python

"""PyMotif

This is a test program.

Description here.

Laurens Bronwasser, lmbronwa@cs.vu.nl
Martijn Vermaat mvermaat@cs.vu.nl

"""

__version__ = "0.1"
__date__ = "2005/10/05"


import sys
from optparse import OptionParser
from Bio import Fasta


def main():

    parser = OptionParser(usage = "usage: %prog -i FILE [options]",
                          version = "PyMotif %s (%s)" % (__version__, __date__),
                          description = "PyMotif reads an input file in Fasta"
                          " format and prints the gene names.")

    parser.add_option("-i", "--input", dest="input", metavar="FILE",
                      help="read FILE in Fasta format")
    parser.add_option("-s", "--samples", dest="samples", action="store_true",
                      default=False, help="print first ten bases of each gene")

    (options, args) = parser.parse_args()

    if not options.input:
        parser.error("input file required")

    # Sample data in Fasta format
    try:
        file = open(options.input)
    except IOError:
        parser.error("could not read file %s" % options.input)

    # Parser for Fasta format
    parser = Fasta.RecordParser()

    # Iterator for sample data
    iterator = Fasta.Iterator(file, parser)

    for record in iterator:
        # Available: record.title, record.sequence, record._colwidth
        print record.title
        if options.samples:
            print "  (First ten bases: %s)" % record.sequence[:10]


if __name__ == "__main__":
    main()
