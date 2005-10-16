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
* The real challenge is not to find a 'good' motif, but to find the best motif
  and to optimize it (convergence algorithm, phase shifting, etc)
* Run the algorithm more than once and show the best result

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
from random import random, randint, choice
from math import log
from optparse import OptionParser
from Bio import Fasta


def main():

    """
    Main program.
    """

    # TODO: put all Gibbs stuff in a gibbs module

    sequences, motif_width, pseudocounts_weight = initialize()

    pseudocounts = calculate_pseudocounts(sequences, pseudocounts_weight)

    random_motif_positions(sequences, motif_width)

    entropies = []
    scores = []

    for i in range(1, 120):

        #print "Step %i:" % i

        # Pick a random sequence
        current_sequence = choice(sequences)
        sequences_minus_current = filter(lambda s: s != current_sequence, sequences)

        #print "Choose sequence with position:", current_sequence['motif_position']

        motif = calculate_motif(sequences_minus_current, motif_width, pseudocounts)

        #print_motif(motif)

        entropies.append(calculate_entropy(motif, motif_width, pseudocounts))
        scores.append(calculate_score(motif, motif_width))

        calculate_position(motif, motif_width, current_sequence)

        #print_sequences(sequences, motif_width)

        # Do a phase shift every 15 iterations (TODO: provide command line option)
        if i % 15 == 0:
            # TODO: provide command line option for max shift (default 0)
            shift = phase_shift(sequences, motif_width, pseudocounts, 2)
            if shift != 0:
                print i, shift

    #print_entropies(entropies)
    #print_scores(scores)


def phase_shift(sequences, motif_width, pseudocounts, shift):

    """
    Calculate motif for some small phase shifts in each sequence and apply
    the fase shift with best entropy.
    """

    # Check if we can really do these shifts
    # TODO: make a difference between shift_left and shift_right
    for s in sequences:
        if s['motif_position'] < shift:
            shift = s['motif_position']
        if len(s['sequence']) - (s['motif_position'] + motif_width) < shift:
            shift = len(s['sequence']) - (s['motif_position'] + motif_width)

    if shift < 1:
        return 0

    # Calculate entropy for no shift
    motif = calculate_motif(sequences, motif_width, pseudocounts)
    old_entropy = calculate_entropy(motif, motif_width, pseudocounts)

    # Start at position -shift
    for s in sequences:
        s['motif_position'] -= shift

    # Keep best shift information
    best_shift = 0
    best_entropy = old_entropy

    # For all shifts in range (-shift, shift)
    for i in range(shift * 2 + 1):

        # We already calculated the no shift case
        if i != 2:

            motif = calculate_motif(sequences, motif_width, pseudocounts)
            entropy = calculate_entropy(motif, motif_width, pseudocounts)

            if entropy < best_entropy:
                best_shift = i - 2
                best_entropy = entropy

        # Shift sequences
        for s in sequences:
            s['motif_position'] += 1

    # Restore to best shift
    # TODO: pick random shift based on distribution
    for s in sequences:
        s['motif_position'] -= (shift - best_shift + 1)

    return best_shift


def calculate_pseudocounts(sequences, weight):

    """
    Return for each base a weighted pseudocount.
    """

    # TODO: use the background frequency for this
    counts = {}
    for base in "ATCG":
        counts[base] = weight * 0.25

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
    Calculate the position weight matrix for the motif for the given
    sequences and their alignments.
    """

    # Populate the matrix with pseudocounts for all positions
    motif = {}
    for base in "ATCG":
        motif[base] = [pseudocounts[base]] * motif_width

    pseudocounts_total = sum(pseudocounts.values())

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
            motif[i][j] /= float(len(sequences)) + pseudocounts_total

    return motif


def calculate_position(motif, motif_width, sequence):

    """
    Calculate new position of motif in sequence based on frequency matrix
    motif with width motif_width.
    """

    # This will be the probability distribution of all positions
    probabilities = []

    # For every word of length motif_width in sequence
    for r in range(len(sequence['sequence']) - motif_width + 1):

        # p_motif: probability that word is generated from motif
        # p_background: probability that word is generated from background
        p_motif = p_background = 1

        # For every position of the motif
        for x in range(motif_width):

            # Multiply by position weight of current base in sequence as
            # defined in the motif matrix
            p_motif *= motif[ sequence['sequence'][r+x] ][x]

            # Multiply by background value of current base in sequence
            # TODO: use background_value[ s[r+x] ]
            p_background *= 0.25

        # Corrected probability for background probability
        probabilities.append(p_motif / p_background)

    # Choose a random position based on probability distribution
    sequence['motif_position'] = choose_random(probabilities)


def calculate_entropy(motif, motif_width, pseudocounts):

    """
    We calculate the relative entropy of the given motif position weight
    matrix, using the background pseudocounts.
    """

    entropy = 0

    for base in "ATCG":
        for i in range(motif_width):
            entropy -= motif[base][i] * log(motif[base][i], 2)

            # I think there'se no need to use the background pseudocounts here,
            # because we already used them in the motif matrix.
            #entropy += motif[base][i] * log(motif[base][i] / pseudocounts[base], 2)

    return entropy


def calculate_score(motif, motif_width):

    score = 0

    for i in range(motif_width):
        mx = 0
        for base in "ATCG":
            if motif[base][i] > mx: mx = motif[base][i]
        score += mx

    return 20 - score


def choose_random(distribution):

    """
    Choose a random integer N with 0 <= N < len(distribution) according to
    distribution of values in distribution.
    """

    # Normalize distribution
    total = sum(distribution)
    distribution = map(lambda n: float(n) / total, distribution)

    # Generate random float between 0 and 1
    position = random()

    current = 0
    for p in range(len(distribution)):
        if distribution[p] + current >= position:
            return p
        current += distribution[p]

    print "Error in choose_random."
    sys.exit(1)


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


def print_entropies(entropies):

    from rpy import r   # TODO: make this import optional

    r.postscript("entropies.ps")
    r.plot(entropies, type='b', xlab="Iterations", ylab="Entropy")
    r.dev_off()


def print_scores(scores):

    # This metric seems to be a little bit less useful than entropy, but it is
    # twice as fast.

    from rpy import r

    r.postscript("scores.ps")
    r.plot(scores, type='b', xlab="Iterations", ylab="Score")
    r.dev_off()


if __name__ == "__main__":
    main()
