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

# print motif occurrence in each sequence
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
* Would it make sense to check for multiple occurrences of the motif in one
  sequence?
* Look at pseudocounts
* Look at gibbs sampling on page 412 of bioinformatics algorithms
* Make a score(motif) function
* The real challenge is not to find a 'good' motif, but to find the best motif
  and to optimize it (convergence algorithm, phase shifting, etc)
* Run the algorithm more than once and show the best result
* In documentation show plots of what happens with phase shifts
  (less stable, but more probable to find better neighbouring alignment)
* heuristic for initial motif positions:
  find random h-word in one sequence, now find the exact same h-word in all
  others and align them
  preferably the h-word should contain highly significant bases (compared to
  background which we already calculated)

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

ITERATIONS_DEFAULT = 80
PHASE_SHIFTS_DEFAULT = 0
PS_FREQUENCY_DEFAULT = 15
PSEUDOCOUNTS_WEIGHT_DEFAULT = 0.1
INIT_NUM_OCCURRENCES_DEFAULT = 0
INIT_PATTERN_WIDTH_DEFAULT = 0


import sys
from gibbs import Gibbs, GibbsError
from optparse import OptionParser
from Bio import Fasta
from random import choice


def main():

    """
    Main program.
    """

    data = initialize()

    g = Gibbs(sequences=data['sequences'],
              motif_width=data['width'],
              pseudocounts_weight=data['weight'])

    g.find_motif(iterations=data['iterations'],
                 phase_shifts=data['shifts'],
                 ps_frequency=data['ps_freq'],
                 initial_num_occurrences=data['init_occurrences'],
                 initial_pattern_width=data['init_width'])

    print_sequences(data['sequences'], data['width'])

    return


def initialize():

    """
    Parses command line options, and reads input Fasta file.

    The returned dictionary contains the following fields:

        sequences   a list of dictionary objects having 'title', 'sequence',
                    and 'motif_position' attributes
        width       width of motif to find
        weight      weight to use for pseudocounts
        iterations  number of non-improving iterations before stopping
        shifts      maximum phase shifts to detect
    """

    # TODO: change description (story about best option settings depending on
    # input genes)
    parser = OptionParser(usage = "usage: %prog -i FILE -w WIDTH [-h] "
                          "[options]",
                          version = "PyMotif %s (%s)" % (VERSION, DATE),
                          description = "PyMotif reads an input file in Fasta"
                          " format and prints the gene names.")

    parser.add_option("-i", "--input", dest="input", metavar="FILE",
                      help="read FILE in Fasta format")
    parser.add_option("-w", "--width", dest="width", metavar="WIDTH",
                      type="int", help="find motif of width WIDTH")
    parser.add_option("-t", "--iterations", dest="iterations",
                      metavar="ITERATIONS", default=ITERATIONS_DEFAULT,
                      type="int", help="number of non-improving iterations ("
                      + str(ITERATIONS_DEFAULT) + ")")
    parser.add_option("-p", "--pseudo", dest="pseudo", metavar="WEIGHT",
                      default=PSEUDOCOUNTS_WEIGHT_DEFAULT, type="float",
                      help="use WEIGHT for weight of pseudocounts (" +
                      str(PSEUDOCOUNTS_WEIGHT_DEFAULT) + ")")
    parser.add_option("-s", "--phase-shifts", dest="shifts", metavar="SHIFTS",
                      default=PHASE_SHIFTS_DEFAULT, type="int",
                      help="detect phase shifts of width SHIFTS (" +
                      str(PHASE_SHIFTS_DEFAULT) + ")")
    parser.add_option("-f", "--ps-frequency", dest="frequency",
                      metavar="FREQ", default=PS_FREQUENCY_DEFAULT,
                      type="int", help="if SHIFTS>0, detect phase shifts "
                      "every FREQ iterations (" + str(PS_FREQUENCY_DEFAULT) +
                      ")")
    parser.add_option("-n", "--init-num-occurrences", dest="initoccurrences",
                      metavar="OCCURRENCES",
                      default=INIT_NUM_OCCURRENCES_DEFAULT, type="int",
                      help="number of base occurrences to use for initial "
                      "positions heuristic (" +
                      str(INIT_NUM_OCCURRENCES_DEFAULT) + ")")
    parser.add_option("-v", "--init-pattern-width", dest="initwidth",
                      metavar="WIDTH", default=INIT_PATTERN_WIDTH_DEFAULT,
                      type="int", help="if OCCURRENCES>0, width of pattern "
                      "to use for initial positions heuristic (same as motif "
                      "width by default)")
    parser.add_option("-c", "--cow", action="store_true", dest="cow",
                      default=False, help="display cow (not recommended)")

    (options, args) = parser.parse_args()

    if options.cow:
        s = ""
        for _ in range(10):
            s += choice("ATCG")
        # Created with the cowsay program
        print """ ____________
< %s >
 ------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\\
                ||----w |
                ||     ||""" % s
        sys.exit(0)

    if not options.input:
        parser.error("input file required")

    if not options.width:
        parser.error("width argument required")

    if options.width < 2:
        parser.error("please use a sane motif width")

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

    # We could do some more error checking on the input file here, like
    # checking there's only ATCG and at least a few of them, but for now
    # this is enough
    if len(sequences) < 2:
        parser.error("found %i sequences in input file %s" % (len(sequences),
                                                              options.input))

    return {'sequences':        sequences,
            'width':            options.width,
            'weight':           options.pseudo,
            'iterations':       options.iterations,
            'shifts':           options.shifts,
            'ps_freq':          options.frequency,
            'init_occurrences': options.initoccurrences,
            'init_width':       options.initwidth}


def print_sequences(sequences, motif_width):

    """
    Print for the occurrence of the motif in each sequence.
    """

    print "Motif occurrences in sequences follow"

    for i in range(len(sequences)):
        start, end = (sequences[i]['motif_position'],
                      sequences[i]['motif_position'] + motif_width)
        print "Sequence #%i  %s  (at position %i)" % (
            i + 1,
            sequences[i]['sequence'][start:end],
            sequences[i]['motif_position'])

    return


if __name__ == "__main__":
    main()
