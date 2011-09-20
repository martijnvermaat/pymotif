#!/usr/bin/env python


"""
PyMotif

Martijn Vermaat <martijn@vermaat.name>

PyMotif is an implementation of the Gibbs sampling algorithm for finding
local alignments of DNA sequences.

Consult the accompanied README file for usage instructions (or run PyMotif
with the -h option) and the documentation directory for implementation
details.

See the file LICENSE for copyright information and the terms and conditions
for copying, distribution and modification of PyMotif.

TODO: readall on repository (in november)
"""


VERSION = "0.1.1"
DATE = "2005/10/25"

ITERATIONS_DEFAULT = 80
PHASE_SHIFTS_DEFAULT = 0
PS_FREQUENCY_DEFAULT = 20
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
    Parse command line options, and read input Fasta file.

    Construct a dictionary contains the following fields:

        sequences         a list of dictionary objects having 'title',
                          'sequence', and 'motif_position' attributes (see
                          also the docstring of gibbs.Gibbs.__init__)
        width             width of motif to find
        weight            weight to use for pseudocounts
        iterations        number of non-improving iterations before stopping
        shifts            maximum phase shifts to detect
        ps_freq           frequency of detecting phase shifts
        init_occurrences  number of base occurrences to use for initial motif
                          positions heuristic
        init_width        width of patterns to use for initial motif positions
                          heuristic

    Return the constructed dictionary.
    """

    parser = OptionParser(usage = "usage: %prog -i FILE -w WIDTH [-h] "
                          "[options]",
                          version = "PyMotif %s (%s)" % (VERSION, DATE),
                          description = "PyMotif is an implementation of the "
                          "Gibbs sampling algorithm for finding local "
                          "alignments of DNA sequences. "
                          "See the accompanied README file for usage "
                          "instructions and the documentation directory for "
                          "implementation details.")

    parser.add_option("-i", "--input", dest="input", metavar="FILE",
                      help="read FILE in Fasta format")
    parser.add_option("-w", "--width", dest="width", metavar="WIDTH",
                      type="int", help="find motif of width WIDTH")
    parser.add_option("-t", "--iterations", dest="iterations",
                      metavar="ITERATIONS", default=ITERATIONS_DEFAULT,
                      type="int", help="number of non-improving iterations "
                      "(default " + str(ITERATIONS_DEFAULT) + ")")
    parser.add_option("-p", "--pseudo", dest="pseudo", metavar="WEIGHT",
                      default=PSEUDOCOUNTS_WEIGHT_DEFAULT, type="float",
                      help="use WEIGHT for weight of pseudocounts (default " +
                      str(PSEUDOCOUNTS_WEIGHT_DEFAULT) + ")")
    parser.add_option("-s", "--phase-shifts", dest="shifts", metavar="SHIFTS",
                      default=PHASE_SHIFTS_DEFAULT, type="int",
                      help="detect phase shifts of width SHIFTS (default " +
                      str(PHASE_SHIFTS_DEFAULT) + ")")
    parser.add_option("-f", "--ps-frequency", dest="frequency",
                      metavar="FREQ", default=PS_FREQUENCY_DEFAULT,
                      type="int", help="if SHIFTS>0, detect phase shifts "
                      "every FREQ iterations (default " +
                      str(PS_FREQUENCY_DEFAULT) + ")")
    parser.add_option("-n", "--init-num-occurrences", dest="initoccurrences",
                      metavar="OCCURRENCES",
                      default=INIT_NUM_OCCURRENCES_DEFAULT, type="int",
                      help="number of base occurrences to use for initial "
                      "positions heuristic (default " +
                      str(INIT_NUM_OCCURRENCES_DEFAULT) + ")")
    parser.add_option("-v", "--init-pattern-width", dest="initwidth",
                      metavar="WIDTH", default=INIT_PATTERN_WIDTH_DEFAULT,
                      type="int", help="if OCCURRENCES>0, width of pattern "
                      "to use for initial positions heuristic (defaults to "
                      "value of --width)")
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
    Print the occurrence of the motif in each sequence.
    """

    print "Motif occurrences in sequences follow"

    for i in range(len(sequences)):
        start, end = (sequences[i]['motif_position'],
                      sequences[i]['motif_position'] + motif_width)
        print "Sequence #%2i  %s  (at position %i)" % (
            i + 1,
            sequences[i]['sequence'][start:end],
            sequences[i]['motif_position'] + 1)

    return


if __name__ == "__main__":
    main()
