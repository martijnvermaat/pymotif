#!/usr/bin/env python


"""
PyMotif

This is a test program.

Description here.


Gibbs sampling algorithm
------------------------

# algorithm from lect4.pdf:
# http://ibivu.cs.vu.nl/teaching/a4g/materials/lect4.pdf

W = width of motif to find (command line argument)
S = set of sequences to analyse (from Fasta file)
b[s] = position of motif in sequence s
M = 4xW matrix of frequencies of bases in motif

# initial (random) values for b
for s in S:
    b[s] = random(0, length(s))

# until convergence calculate new M and b
while True:

    cur_s = random s in S

    M = calculate_m(W, S - cur_s, b)
    b[cur_s] = calculate_position(M, cur_s)

    # how to detect convergence?
    if converged: break

# now we have 'good' M and b

# print motif occurence in each sequence
for s in S:
    print s[ b[s]:(W-b[s] ]

# M can be used to highlight matching bases (as done by AlignAce)

def calculate_m(W, S, b):

    # S = set of sequences
    # n = number of sequences
    # b[s] = position of motif in sequence s
    # M = 4xW matrix of frequencies of bases in motif

    n = length(S)

    # for all bases and all positions of motif
    for i in ['A','T','C','G']:
        for j in [0..W]:

            # M[i][j] = frequency of i on position j
            M[i][j] = 0
            for s = S:
                if s[ b[s]+j-1 ] == i: M[i][j]++
                M[i][j] *= 1/n

    # implementation from page 6 of lect4.pdf, slide 4

    # this implementation is very naive, because we shouldn't
    # calculate entire M every iteration. instead do something
    # with 'pseudocounts'.
    # todo: what are pseudocounts?

    return M

def calculate_position(M, s):

    # M = 4xW matrix of frequencies of bases in motif
    # W = width of motif
    # s = sequence to find motif in
    # p = position of motif in s

    # todo: look at page 7 of lect4.pdf, slide 2

    return p


Laurens Bronwasser, lmbronwa@cs.vu.nl
Martijn Vermaat mvermaat@cs.vu.nl
"""


VERSION = "0.1"
DATE = "2005/10/05"


import sys
from optparse import OptionParser
from Bio import Fasta


def main():

    parser = OptionParser(usage = "usage: %prog -i FILE [options]",
                          version = "PyMotif %s (%s)" % (VERSION, DATE),
                          description = "PyMotif reads an input file in Fasta"
                          " format and prints the gene names.")

    parser.add_option("-i", "--input", dest="input", metavar="FILE",
                      help="read FILE in Fasta format")
    parser.add_option("-s", "--samples", dest="samples", action="store_true",
                      default=False, help="print first ten bases of each gene")
    parser.add_option("-w", "--width", dest="width", metavar="WIDTH",
                      help="find motif of width WIDTH")

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
