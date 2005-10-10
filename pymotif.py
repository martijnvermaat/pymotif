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
    b[s] = random(0, length(s) - W)

# until convergence calculate new M and b
while True:

    cur_s = random s in S
    # alt: cur_s = round robin s in S
    # we could do a performance check on this alternative
    # not important at the moment

    M = calculate_m(W, S - cur_s, b)
    b[cur_s] = calculate_position(W, M, cur_s)

    # how to detect convergence?
    if converged: break

# now we have 'good' M and b

# print motif occurence in each sequence
for s in S:
    print s[ b[s]:(W-b[s]) ]

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
            for s in S:
                if s[ b[s]+j-1 ] == i: M[i][j]++
            M[i][j] *= 1/n

    # implementation from page 6 of lect4.pdf, slide 4

    # this implementation is very naive, because we shouldn't
    # calculate entire M every iteration. instead do something
    # with 'pseudocounts'.
    # todo: what are pseudocounts?

    return M

def calculate_position(W, M, s):

    # M = 4xW matrix of frequencies of bases in motif
    # W = width of motif
    # s = sequence to find motif in
    # r = current position of motif in s to consider
    # p = ultimate position of motif in s
    # C = list of probabilities for each position of motif in s

    # todo: look at page 7 of lect4.pdf, slide 2

    # i don't completely understand this slide
    # * what do the Qr and Pr say?
    # * in the computation of Pr, where do the p's come from?
    #   (later: i think little p is the same as our M)
    #   (later: i don't think so, p's have something to do with background)
    # * i guess for every K letter word in s, the r index is the position
    #   of this word in s?
    # * by K she actually means W
    # * great chance with K she means W :)

    # ok, according to my mail at 00:45, 10 oct I suggest the following
    # implementation (comments are left out, consult the mail)

    # for every word of length W in s
    for r in range(len(s) - W):

        Qr = Pr = 1

        # for every position of the motif
        for x in range(W):

            # multiply by position weight of current base in s as defined in
            # the matrix
            Qr *= M[ s[r+x] ][x]

            # multiply by background value of current base in s
            Pr *= background_value[ s[r+x] ]

        # corrected Qr value
        # devide by zero?
        C[r] = Qr / Pr

    # pick random position according to distribution of C[r] values
    p = r for maximum C[r]

    return p


Random notes:
* Would it make sense to check for multiple occurences of the motif in one
  sequence?
* Look at pseudocounts
* Look at gibbs sampling on page 412 of bioinformatics algorithms


Laurens Bronwasser, lmbronwa@cs.vu.nl
Martijn Vermaat mvermaat@cs.vu.nl
"""


VERSION = "0.1"
DATE = "2005/10/10"


import sys
from random import randint
from optparse import OptionParser
from Bio import Fasta


def main():

    """
    Main program.

    Real algorithm should be in a separate gibbs() function.
    """

    sequences, motif_width = initialize()

    random_motif_positions(sequences, motif_width)

    motif = calculate_motif(sequences, motif_width)

    # Print start of each sequence and current motif position
    for s in sequences:
        print "%s... motif at position %s" % (s['sequence'][:40], s['motif_position'])

    # Print position weight matrix for motif
    for base in "ATCG":
        print base,
        for weight in motif[base]:
            #print "%1.2f" % weight,   # we have floats
            print "%3i" % weight,      # we have ints
        print


def random_motif_positions(sequences, motif_width):

    """
    Populate the list of sequences with a random position of the motif for
    each sequence.
    """

    for s in sequences:
        s['motif_position'] = randint(0, len(s['sequence']) - motif_width)


def calculate_motif(sequences, motif_width):

    """
    Calculate the position weight matrix for the motif.

    This implementation is very naive, use pseudocounts.
    """

    # Populate the matrix with 0 for all positions
    motif = {'A': [0 for _ in range(motif_width)],
             'T': [0 for _ in range(motif_width)],
             'C': [0 for _ in range(motif_width)],
             'G': [0 for _ in range(motif_width)]}

    # For all bases and all positions of motif
    for i in "ATCG":
        for j in range(motif_width):

            # For each sequence, add 1 if it has base i at position j
            for s in sequences:
                position = s['motif_position'] + j
                if s['sequence'][position] == i:
                    motif[i][j] += 1

            # Divide by the number of sequences and we have the weight of base
            # i at position j

            # I think we shouldn't do this devision
            # Perhaps only on output of the matrix to the screen
            # (Have to check if this is ok for rest of algorithm)
            #motif[i][j] /= float(len(sequences))

    return motif


def initialize():

    """
    Parse command line options, and read input Fasta file.

    Provides:
      sequences: a list of dictionary objects having 'title', 'sequence', and
                 'motif_position' attributes.
      width: width of motif to find

    Return these as a pair (sequences, width).
    """

    parser = OptionParser(usage = "usage: %prog -i FILE -w  WIDTH [options]",
                          version = "PyMotif %s (%s)" % (VERSION, DATE),
                          description = "PyMotif reads an input file in Fasta"
                          " format and prints the gene names.")

    parser.add_option("-i", "--input", dest="input", metavar="FILE",
                      help="read FILE in Fasta format")
    parser.add_option("-w", "--width", dest="width", metavar="WIDTH",
                      type="int", help="find motif of width WIDTH")
    parser.add_option("-f", "--format", action="store_true", dest="format",
                      default=False, help="format teacher's harddrive")

    (options, args) = parser.parse_args()

    if not options.input:
        parser.error("input file required")

    if not options.width:
        parser.error("width argument required")

    # Read contents of Fasta file
    try:
        file = open(options.input)
    except IOError:
        parser.error("could not read file %s" % options.input)

    fasta_parser = Fasta.RecordParser()

    # Iterator for sample data
    fasta_iterator = Fasta.Iterator(file, fasta_parser)

    # A list containing a dictionary object for each sequence
    sequences = [{'title':          record.title,
                  'sequence':       record.sequence,
                  'motif_position': 0}
                 for record in fasta_iterator]

    return sequences, options.width


if __name__ == "__main__":
    main()
