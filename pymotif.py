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
* Make a score(motif) function

Pointers:
http://ibivu.cs.vu.nl/teaching/a4g/materials/lect3-handout.pdf
http://ibivu.cs.vu.nl/teaching/a4g/materials/lect4.pdf
http://ai.stanford.edu/~serafim/CS262_2005/LectureNotes/Lecture15.pdf
http://www.snubi.org/~jhohn2/20020611.ppt
http://genomics10.bu.edu/ulask/mcmc/mcmc_cs591.pdf
http://melina.hgc.jp/expl_par.html


Martijn Vermaat mvermaat@cs.vu.nl
"""


VERSION = "0.1"
DATE = "2005/10/12"

PSEUDOCOUNTS_DEFAULT_WEIGHT = 0.1


import sys
from random import randint, choice
from optparse import OptionParser
from Bio import Fasta


def main():

    """
    Main program.

    Real algorithm should be in a separate gibbs() function.
    """

    sequences, motif_width, pseudocounts_weight = initialize()

    pseudocounts = calculate_pseudocounts(sequences, pseudocounts_weight)

    random_motif_positions(sequences, motif_width)

    print "Initialization:"

    motif = calculate_motif(sequences, motif_width, pseudocounts)

    print_motif(motif)

    print_sequences(sequences, motif_width)

    for i in range(200):

        print "Step %i:" % i

        # Pick a random sequence
        current_sequence = choice(sequences)
        sequences_minus_current = filter(lambda s: s != current_sequence, sequences)

        print "Choose sequence with position:", current_sequence['motif_position']

        motif = calculate_motif(sequences_minus_current, motif_width, pseudocounts)

        print_motif(motif)

        calculate_position(motif, motif_width, current_sequence)

        print_sequences(sequences, motif_width)


def calculate_pseudocounts(sequences, weight):

    """
    Return for each base a weighted pseudocount.
    """

    # TODO: use the background frequency for this
    counts = {'A': weight * 0.25,
              'T': weight * 0.25,
              'C': weight * 0.25,
              'G': weight * 0.25}

    return counts

def random_motif_positions(sequences, motif_width):

    """
    Populate the list of sequences with a random position of the motif for
    each sequence.
    """

    for s in sequences:
        s['motif_position'] = randint(0, len(s['sequence']) - motif_width)


def calculate_motif(sequences, motif_width, pseudocounts):

    """
    Calculate the position weight matrix for the motif.
    """

    # TODO: this implementation is very naive

    # Populate the matrix with pseudocounts for all positions
    motif = {'A': [pseudocounts['A']] * motif_width,
             'T': [pseudocounts['T']] * motif_width,
             'C': [pseudocounts['C']] * motif_width,
             'G': [pseudocounts['G']] * motif_width}

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


def calculate_position(motif, motif_width, sequence):

    """
    Calculate new position of motif in sequence.
    """

    # TODO: pick random based on probabilities instead of the highest
    highest_probability = 0
    position = 0

    # for every word of length W in s
    for r in range(len(sequence['sequence']) - motif_width + 1):

        Qr = Pr = 1

        # for every position of the motif
        for x in range(motif_width):

            # multiply by position weight of current base in sequence as
            # defined in the motif matrix
            Qr *= motif[ sequence['sequence'][r+x] ][x]

            # multiply by background value of current base in sequence
            Pr *= 0.25 # TODO: use background_value[ s[r+x] ]

        # corrected Qr value
        # devide by zero?
        #sample_probabilities.append(Qr / Pr)
        probability = Qr / Pr
        if probability > highest_probability:
            highest_probability = probability
            position = r

    sequence['motif_position'] = position

    print "New position:", position

    return position


def initialize():

    """
    Parse command line options, and read input Fasta file.

    Provides:
      sequences: a list of dictionary objects having 'title', 'sequence', and
                 'motif_position' attributes.
      width: width of motif to find
      weight: weight of pseudocounts

    Return these as a tuple (sequences, width, weight).
    """

    # TODO: change description
    parser = OptionParser(usage = "usage: %prog -i FILE -w WIDTH [options]",
                          version = "PyMotif %s (%s)" % (VERSION, DATE),
                          description = "PyMotif reads an input file in Fasta"
                          " format and prints the gene names.")

    parser.add_option("-i", "--input", dest="input", metavar="FILE",
                      help="read FILE in Fasta format")
    parser.add_option("-w", "--width", dest="width", metavar="WIDTH",
                      type="int", help="find motif of width WIDTH")
    parser.add_option("-p", "--pseudo", dest="pseudo", metavar="WEIGHT",
                      default=PSEUDOCOUNTS_DEFAULT_WEIGHT, type="float",
                      help="use WEIGHT for weight of pseudocounts (defaults to " +
                      str(PSEUDOCOUNTS_DEFAULT_WEIGHT) + ")")
    parser.add_option("-f", "--format", action="store_true", dest="format",
                      default=False, help="format teacher's harddrive (not recommended)")

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

    return sequences, options.width, options.pseudo


def print_motif(motif):
    # Print position weight matrix for motif
    for base in "ATCG":
        print base,
        # Do not uncomment this comment
        for weight in motif[base]:
            print "%1.2f" % weight,   # we have floats
        print


def print_sequences(sequences, motif_width):
    # Print motif occurence in each sequence and motif position
    for s in sequences:
        start, end = s['motif_position'], s['motif_position'] + motif_width
        print "%s... motif at position %s" % (s['sequence'][start:end], s['motif_position'])


if __name__ == "__main__":
    main()
