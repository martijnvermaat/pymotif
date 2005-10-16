"""
Set of functions implementing the Gibbs sampling algorithm to find local
alignments of sequences.


TODO: use __ in private names


Martijn Vermaat mvermaat@cs.vu.nl
"""


import random
import math


class Gibbs:


    """
    This class implements the Gibbs sampling algorithm.
    """


    __sequences = []
    __motif_width = 0
    __pseudocounts = {}


    def __init__(self, sequences, motif_width, pseudocounts_weight):

        """
        Gibbs takes a list of sequence objects, the width of the motif to find
        and the weight to use for pseudocounts.

        The sequence object are dictionaries with keys:

            title           title of the sequence (e.g. gene name)
            sequence        string with ATCG occurences in sequence
            motif_position  integer position of the motif in sequence
        """

        self.__sequences = sequences
        self.__motif_width = motif_width

        self.__calculate_pseudocounts(pseudocounts_weight)

        return


    def find_motif(self, runs=1, iterations=50, phase_shifts=0):

        """
        Finds a motif in sequences using Gibbs sampling.

        Parameters:

            runs          times to run algorithm
            iterations    non-improving iterations before quiting
            phase_shifts  number of phase shifts to examine

        Stores the best alignment found.
        """

        overall_best_entropy = 9999
        overall_best_alignment = []

        for r in range(runs):

            # TODO: maybe we could make some heuristics up for the next
            # initialization instead of random positions every run
            self.__random_motif_positions()

            best_entropy = 9999
            best_alignment = []

            i = 0

            while i < iterations:

                i += 1

                # Pick a random sequence
                current_sequence = random.choice(self.__sequences)
                sequences_minus_current = filter(
                    lambda s: s != current_sequence,
                    self.__sequences)

                motif = self.__calculate_motif(sequences_minus_current)

                #print_motif(motif)

                #entropies.append(g.calculate_entropy(motif,
                #                                     motif_width,
                #                                     pseudocounts))

                self.__calculate_position(motif, current_sequence)

                # Do a phase shift every 15 iterations
                # TODO: provide command line option
                if i > 0 and i % 15 == 0:
                    shift = self.__phase_shift(phase_shifts)
                    #if shift != 0:
                        #print "Shift:", i, shift

                entropy = self.__calculate_entropy(motif)
                if entropy < best_entropy:
                    best_entropy = entropy
                    best_alignment = [s['motif_position']
                                      for s in self.__sequences]
                    i = 0

            # Restore to best entropy
            for i in range(len(self.__sequences)):
                self.__sequences[i]['motif_position'] = best_alignment[i]

            motif = self.__calculate_motif(self.__sequences)
            entropy = self.__calculate_entropy(motif)

            #print entropy

            if entropy < overall_best_entropy:
                overall_best_entropy = entropy
                overall_best_alignment = [s['motif_position']
                                          for s in self.__sequences]

        # Restore to best entropy
        for i in range(len(self.__sequences)):
            self.__sequences[i]['motif_position'] = overall_best_alignment[i]

        motif = self.__calculate_motif(self.__sequences)
        print self.__calculate_entropy(motif)

        return


    def __phase_shift(self, shift):

        """
        Calculate motif for some small phase shifts in each sequence and apply
        the fase shift with best entropy.
        """

        # Check if we can really do these shifts
        # TODO: make a difference between shift_left and shift_right
        for s in self.__sequences:
            if s['motif_position'] < shift:
                shift = s['motif_position']
            right = len(s['sequence']) - (s['motif_position']
                                          + self.__motif_width)
            if  right < shift:
                shift = right

        if shift < 1:
            return 0

        # Calculate entropy for no shift
        motif = self.__calculate_motif(self.__sequences)
        old_entropy = self.__calculate_entropy(motif)

        # Start at position -shift
        for s in self.__sequences:
            s['motif_position'] -= shift

        # Keep best shift information
        best_shift = 0
        best_entropy = old_entropy

        # For all shifts in range (-shift, shift)
        for i in range(shift * 2 + 1):

            # We already calculated the no shift case
            if i != 2:

                motif = self.__calculate_motif(self.__sequences)
                entropy = self.__calculate_entropy(motif)

                if entropy < best_entropy:
                    best_shift = i - 2
                    best_entropy = entropy

            # Shift sequences
            for s in self.__sequences:
                s['motif_position'] += 1

        # Restore to best shift
        # TODO: pick random shift based on distribution
        for s in self.__sequences:
            s['motif_position'] -= (shift - best_shift + 1)

        return best_shift


    def __calculate_pseudocounts(self, weight):

        """
        Calculate for each base a weighted pseudocount.
        """

        # TODO: use the background frequency for this
        for base in "ATCG":
            self.__pseudocounts[base] = weight * 0.25

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


    def __calculate_motif(self, sequences):

        """
        Calculate the position weight matrix of the motif for the given
        sequences and their alignments.

        The resulting matrix is a dictionary with for each key in ATCG a list
        of length motif_width with position weights.
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
                # TODO: use background_value[ s[r+x] ]
                p_background *= 0.25

            # Corrected probability for background probability
            probabilities.append(p_motif / p_background)

        # Choose a random position based on probability distribution
        sequence['motif_position'] = self.__choose_random(probabilities)

        return


    def __calculate_entropy(self, motif):

        """
        Calculate the relative entropy of the given motif matrix.
        """

        entropy = 0

        for base in "ATCG":
            for i in range(self.__motif_width):
                entropy -= motif[base][i] * math.log(motif[base][i], 2)

                # I think there'se no need to use the background pseudocounts
                # here, because we already used them in the motif matrix.
                #entropy += motif[base][i] * log(motif[base][i] /
                #                                self.__pseudocounts[base], 2)

        return entropy


    def __choose_random(self, distribution):

        """
        Choose a random integer N with 0 <= N < len(distribution) according to
        distribution of values in distribution.
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


class GibbsError(Exception):
    """
    Exception thrown by Gibbs instances.
    """
    # Our own exception doesn't need anything fancy.
    pass
