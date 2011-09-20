"""
This module contains an implementation of the Gibbs sampling algorithm used
for finding local alignments of DNA sequences.

Classes:

    Gibbs       Implementation of Gibbs sampling
    GibbsError  Class of errors raised by the Gibbs sampling algorithm

Martijn Vermaat <martijn@vermaat.name>

See the file LICENSE for copyright information and the terms and conditions
for copying, distribution and modification of PyMotif.
"""


import random
import math


PS_ITERATIONS_FACTOR = 4
ONE_PLUS_ONE = 2

DEBUG = False
DEBUG_ENTROPIES_FILE = "entropies.ps"


class Gibbs:


    """
    An implementation of the Gibbs sampling algorithm.

    Methods:

        find_motif  start algorithm on provided DNA sequences
    """


    __sequences = []
    __motif_width = 0
    __pseudocounts = {}


    def __init__(self, sequences, motif_width, pseudocounts_weight=0.1):

        """
        Create a new Gibbs instance with provided sequences and motif width.
        Call find_motif after this to run the Gibbs sampling algorithm.

        Parameters:

            sequences            list of sequence objects (described below)
            motif_width          width of the motif to find
            pseudocounts_weight  weight to use for pseudocounts (default 0.1)

        A sequence object is a dictionary with the following keys:

            title           title of the sequence (e.g. gene name)
            sequence        string with ATCG occurrences in sequence
            motif_position  integer position of the motif in sequence
        """

        self.__sequences = sequences
        self.__motif_width = motif_width

        self.__calculate_pseudocounts(pseudocounts_weight)

        return


    def find_motif(self, iterations=50, phase_shifts=0, ps_frequency=12,
                   initial_num_occurrences=0, initial_pattern_width=0):

        """
        Find a motif in sequences using Gibbs sampling.

        Parameters:

            iterations               non-improving iterations before quiting
                                     (default 50)
            phase_shifts             number of phase shifts to examine
                                     (default 0)
            ps_frequency             frequency of iterations detecting phase
                                     shifts (default 12)
            initial_num_occurrences  number of base occurrences in pattern for
                                     heuristic of initial motif positions
                                     (default 0)
            initial_pattern_width    width of pattern for heuristic of initial
                                     motif positions (default value of
                                     motif_width)

        Store the best alignment found in the sequence objects.
        """

        if initial_num_occurrences > 0:

            if initial_pattern_width == 0:
                initial_pattern_width = self.__motif_width

            self.__initial_motif_positions(initial_num_occurrences,
                                           initial_pattern_width)

        else:

            self.__random_motif_positions()

        best_entropy = 0
        best_alignment = [0] * len(self.__sequences)

        if DEBUG: entropies = []

        i = j = 0

        while i < iterations:

            i += 1
            j += 1

            # Pick a random sequence
            current_sequence = random.choice(self.__sequences)
            sequences_minus_current = filter(
                lambda s: s != current_sequence,
                self.__sequences)

            motif = self.__calculate_motif(sequences_minus_current)

            if DEBUG:
                entropies.append(self.__calculate_entropy(motif))

            self.__calculate_position(motif, current_sequence)

            # Do a phase shift every ps_frequency iterations
            if phase_shifts > 0 and j % ps_frequency == 0:
                shift = self.__phase_shift(phase_shifts)
                # Give the algorithm a chance of recovering from
                # badly chosen phase shifts
                # The correction for i should never be bigger than
                # ps_frequency because it would loop forever
                i -= ps_frequency / PS_ITERATIONS_FACTOR

            # Okay, this is based on N-1 sequences, but practice shows this
            # is good enough to remember 'best' alignments
            entropy = self.__calculate_entropy(motif)
            if entropy > best_entropy:
                best_entropy = entropy
                best_alignment = [s['motif_position']
                                  for s in self.__sequences]
                i = 0

        # Restore to best entropy
        for i in range(len(self.__sequences)):
            self.__sequences[i]['motif_position'] = best_alignment[i]

        if DEBUG:
            motif = self.__calculate_motif(self.__sequences)
            entropy = self.__calculate_entropy(motif)
            self.__print_entropies(entropies)
            print "Entropy:", best_entropy
            self.__print_motif(motif)

        return


    def __phase_shift(self, shift):

        """
        Calculate motif for some small phase shifts in each sequence and apply
        the fase shift with best entropy.

        See the PyMotif documentation for more information.

        Parameters:

            shift  maximum width of shifts to consider

        Return the width of the phase shift applied.
        """

        shift_left = shift_right = shift

        # Check if we can really do these shifts
        for s in self.__sequences:
            if s['motif_position'] < shift_left:
                shift_left = s['motif_position']
            right = len(s['sequence']) - (s['motif_position']
                                          + self.__motif_width)
            if  right < shift_right:
                shift_right = right

        if shift_left < 1 and shift_right < 1:
            return 0

        # Calculate entropy for no shift
        motif = self.__calculate_motif(self.__sequences)
        old_entropy = self.__calculate_entropy(motif)

        # Start at left-most position
        for s in self.__sequences:
            s['motif_position'] -= shift_left

        # Keep entropy for all shifts
        entropies = []

        # For all shifts
        for i in range(shift_left + 1 + shift_right):

            # We already calculated the no shift case
            if i == shift_left:
                entropies.append(old_entropy)

            else:
                motif = self.__calculate_motif(self.__sequences)
                entropy = self.__calculate_entropy(motif)
                entropies.append(entropy)

            # Shift sequences
            for s in self.__sequences:
                s['motif_position'] += 1


        # Make distribution based on entropies
        lowest = min(entropies)
        distribution = map(lambda e: e - lowest, entropies)

        best = self.__choose_random(distribution)

        # Restore to best shift
        for s in self.__sequences:
            s['motif_position'] -= (shift_left + 1 + shift_right) - best

        return best - shift_left


    def __calculate_pseudocounts(self, weight):

        """
        Calculate for each base a weighted pseudocount. Store results in
        __pseudocounts dictionary.

        Parameters:

            weight  weight to use for pseudocounts
        """

        # Total number of bases in sequences
        total = float(sum([len(s['sequence']) for s in self.__sequences]))

        for base in "ATCG":
            # Number of occurrences of this base in sequences
            n = sum([s['sequence'].count(base) for s in self.__sequences])
            self.__pseudocounts[base] = weight * (n / total)

        return


    def __random_motif_positions(self):

        """
        Populate the list of sequences with a random position of the motif for
        each sequence.
        """

        for s in self.__sequences:
            s['motif_position'] = random.randint(
                0, len(s['sequence']) - self.__motif_width)

        return


    def __initial_motif_positions(self, number_of_occurrences, pattern_width):

        """
        Populate the list of sequences with a random position of the motif for
        each sequence based on the following heuristic.

        Parameters:

            number_of_occurrences  number of occurrences of base to look for
                                   in patterns
            pattern_width          width of patterns to look for

        The idea is that significant motifs are mostly the ones that contain
        relatively a lot of bases with a low background frequency.
        In order to guide the algorithm to regions in the sequences where more
        of these bases occur, we do the following:

            Pick the base with the lowest background frequency. Now, in each
            sequence, find all words of length pattern_width with at least
            number_of_occurrences occurrences of that base. Choose a random
            one of these position and initialize the motif position of the
            sequence by aligning it with the random chosen position.

        Of course, significant motifs might not always contain a certain
        number of occurrences of the base with lowest background frequency,
        let alone within a certain amount of bases apart from each other.
        """

        # Get base with lowest frequency
        lowest_freq = self.__pseudocounts['A']
        lowest_base = 'A'
        for base in "ATCG":
            if self.__pseudocounts[base] < lowest_freq:
                lowest_freq = self.__pseudocounts[base]
                lowest_base = base

        for s in self.__sequences:

            positions = []

            occurrences = number_of_occurrences

            while occurrences > 0:

                # For every word of length pattern_width in sequence
                for r in range(len(s['sequence']) - pattern_width + 1):

                    # Count occurrences of lowest_base in word
                    n = s['sequence'][r : r+pattern_width].count(lowest_base)

                    if n >= occurrences:
                        positions.append(r)

                if len(positions) > 0:
                    break

                # Try the same with one occurrence of lowest_base less
                occurrences -= 1

            if len(positions) > 0:

                position = random.choice(positions)
                position += pattern_width / 2
                position -= self.__motif_width / 2
                if position < 0:
                    position = 0
                if position > len(s['sequence']) - self.__motif_width:
                    position = len(s['sequence']) - self.__motif_width

            else:

                # This will not happen a lot, but choose a random position if
                # the sequence doesn't have at least one occurrence of
                # lowest_base
                position = random.randint(
                    0, len(s['sequence']) - self.__motif_width)

            s['motif_position'] = position

        return


    def __calculate_motif(self, sequences):

        """
        Calculate the position weight matrix of the motif for the given
        sequences and their alignments.

        Parameters:

            sequences  list of sequence objects to calculate motif for

        Return the matrix as a dictionary with for each key in ATCG a list of
        length motif_width with position weights.
        """

        # Populate the matrix with pseudocounts for all positions
        motif = {}
        for base in "ATCG":
            motif[base] = [self.__pseudocounts[base]] * self.__motif_width

        pseudocounts_total = sum(self.__pseudocounts.values())

        # For all bases and all positions of motif
        for i in "ATCG":
            for j in range(self.__motif_width):

                # For each sequence, add 1 if it has base i at position j
                for s in sequences:
                    position = s['motif_position'] + j
                    if s['sequence'][position] == i:
                        motif[i][j] += 1

                # Divide by the number of sequences and we have the weight of
                # base i at position j
                motif[i][j] /= float(len(sequences)) + pseudocounts_total

        return motif


    def __calculate_position(self, motif, sequence):

        """
        Calculate new position of motif in sequence based on motif matrix.

        Parameters:

            motif     motif to use for calculating new position
            sequence  sequence to align using motif
        """

        # This will be the probability distribution of all positions
        probabilities = []

        # For every word of length motif_width in sequence
        for r in range(len(sequence['sequence']) - self.__motif_width + 1):

            # p_motif: probability that word is generated from motif
            # p_background: probability that word is generated from background
            p_motif = p_background = 1

            # For every position of the motif
            for x in range(self.__motif_width):

                # Multiply by position weight of current base in sequence as
                # defined in the motif matrix
                p_motif *= motif[ sequence['sequence'][r+x] ][x]

                # Multiply by background value of current base in sequence
                p_background *= self.__pseudocounts[
                    sequence['sequence'][r+x]]

            # Corrected probability for background probability
            probabilities.append(p_motif / p_background)

        # Choose a random position based on probability distribution
        sequence['motif_position'] = self.__choose_random(probabilities)

        return


    def __calculate_entropy(self, motif):

        """
        Calculate the relative entropy of the given motif matrix.

        Parameters:

            motif  calculate relative entropy for this motif

        Return the calculated entropy.
        """

        entropy = 0

        for base in "ATCG":
            for i in range(self.__motif_width):
                entropy += motif[base][i] * math.log(
                    motif[base][i] / self.__pseudocounts[base], 2)

        return entropy


    def __choose_random(self, distribution):

        """
        Choose a random integer N with 0 <= N < len(distribution) according to
        distribution of values in distribution.

        Parameters:

            distribution  list of values to use as a distribution

        Return randomly chosen integer.
        """

        # Normalize distribution
        total = sum(distribution)
        distribution = map(lambda n: float(n) / total, distribution)

        # Generate random float between 0 and 1
        position = random.random()

        current = 0
        for p in range(len(distribution)):
            if distribution[p] + current >= position:
                return p
            current += distribution[p]

        # Hmm we should have already left
        raise GibbsError("Error in biased random function (choose_random)")


    def __print_entropies(self, entropies):

        """
        This is for debugging purposes.
        """

        try:
            from rpy import r
        except:
            print "Could not import rpy module"
            return

        r.postscript(DEBUG_ENTROPIES_FILE)
        r.plot(entropies, type='b', xlab="Iterations", ylab="Entropy")
        r.dev_off()

        return


    def __print_motif(self, motif):

        """
        This is for debugging purposes.
        """

        # Print position weight matrix for motif
        for base in "ATCG":
            print base,
            for weight in motif[base]:
                print "%1.2f" % weight,
            print

        # Print significant motif positions
        print '',
        for i in range(self.__motif_width):
            if max([motif[base][i] for base in "ATCG"]) > .6:
                print "   *",
            else:
                print "    ",
        print

        return


    def __print_sequences(self):

        """
        This is for debugging purposes.
        """

        # Print motif occurrence in each sequence and motif position
        for s in self.__sequences:
            start, end = s['motif_position'], s['motif_position'] + self.__motif_width
            print "%s motif at position %s" % (s['sequence'][start:end], s['motif_position'])


class GibbsError(Exception):
    """
    Exception thrown by Gibbs instances.
    """
    # Our own exception doesn't need anything fancy.
    pass
