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


    sequences = []
    motif_width = 0
    pseudocounts = {}


    def __init__(self, sequences, motif_width, pseudocounts_weight):

        """
        Gibbs takes a list of sequence objects, the width of the motif to find
        and the weight to use for pseudocounts.

        The sequence object are dictionaries with keys:

            title           title of the sequence (e.g. gene name)
            sequence        string with ATCG occurences in sequence
            motif_position  integer position of the motif in sequence
        """

        self.sequences = sequences
        self.motif_width = motif_width

        self.calculate_pseudocounts(pseudocounts_weight)

        return


    def find_motif(self, iterations=100):

        """
        Find a motif in sequences applying Gibbs sampling with the given
        number of iterations.
        """

        self.random_motif_positions()

        for i in range(iterations):

            # Pick a random sequence
            current_sequence = random.choice(self.sequences)
            sequences_minus_current = filter(lambda s: s != current_sequence, self.sequences)

            motif = self.calculate_motif(sequences_minus_current)

            #print_motif(motif)

            #entropies.append(g.calculate_entropy(motif, motif_width, pseudocounts))

            self.calculate_position(motif, current_sequence)

            # Do a phase shift every 15 iterations (TODO: provide command line option)
            if i > 0 and i % 15 == 0:
                # TODO: provide command line option for max shift (default 0)
                shift = self.phase_shift(2)
                if shift != 0:
                    print i, shift

        return self.calculate_motif(self.sequences)


    def phase_shift(self, shift):

        """
        Calculate motif for some small phase shifts in each sequence and apply
        the fase shift with best entropy.
        """

        # Check if we can really do these shifts
        # TODO: make a difference between shift_left and shift_right
        for s in self.sequences:
            if s['motif_position'] < shift:
                shift = s['motif_position']
            right = len(s['sequence']) - (s['motif_position'] + self.motif_width)
            if  right < shift:
                shift = right

        if shift < 1:
            return 0

        # Calculate entropy for no shift
        motif = self.calculate_motif(self.sequences)
        old_entropy = self.calculate_entropy(motif)

        # Start at position -shift
        for s in self.sequences:
            s['motif_position'] -= shift

        # Keep best shift information
        best_shift = 0
        best_entropy = old_entropy

        # For all shifts in range (-shift, shift)
        for i in range(shift * 2 + 1):

            # We already calculated the no shift case
            if i != 2:

                motif = self.calculate_motif(self.sequences)
                entropy = self.calculate_entropy(motif)

                if entropy < best_entropy:
                    best_shift = i - 2
                    best_entropy = entropy

            # Shift sequences
            for s in self.sequences:
                s['motif_position'] += 1

        # Restore to best shift
        # TODO: pick random shift based on distribution
        for s in self.sequences:
            s['motif_position'] -= (shift - best_shift + 1)

        return best_shift


    def calculate_pseudocounts(self, weight):

        """
        Calculate for each base a weighted pseudocount.
        """

        # TODO: use the background frequency for this
        for base in "ATCG":
            self.pseudocounts[base] = weight * 0.25

        return


    def random_motif_positions(self):

        """
        Populate the list of sequences with a random position of the motif for
        each sequence.
        """

        for s in self.sequences:
            s['motif_position'] = random.randint(
                0, len(s['sequence']) - self.motif_width)

        return


    def calculate_motif(self, sequences):

        """
        Calculate the position weight matrix of the motif for the given
        sequences and their alignments.

        The resulting matrix is a dictionary with for each key in ATCG a list
        of length motif_width with position weights.
        """

        # Populate the matrix with pseudocounts for all positions
        motif = {}
        for base in "ATCG":
            motif[base] = [self.pseudocounts[base]] * self.motif_width

        pseudocounts_total = sum(self.pseudocounts.values())

        # For all bases and all positions of motif
        for i in "ATCG":
            for j in range(self.motif_width):

                # For each sequence, add 1 if it has base i at position j
                for s in sequences:
                    position = s['motif_position'] + j
                    if s['sequence'][position] == i:
                        motif[i][j] += 1

                # Divide by the number of sequences and we have the weight of
                # base i at position j
                motif[i][j] /= float(len(sequences)) + pseudocounts_total

        return motif


    def calculate_position(self, motif, sequence):

        """
        Calculate new position of motif in sequence based on motif matrix.
        """

        # This will be the probability distribution of all positions
        probabilities = []

        # For every word of length motif_width in sequence
        for r in range(len(sequence['sequence']) - self.motif_width + 1):

            # p_motif: probability that word is generated from motif
            # p_background: probability that word is generated from background
            p_motif = p_background = 1

            # For every position of the motif
            for x in range(self.motif_width):

                # Multiply by position weight of current base in sequence as
                # defined in the motif matrix
                p_motif *= motif[ sequence['sequence'][r+x] ][x]

                # Multiply by background value of current base in sequence
                # TODO: use background_value[ s[r+x] ]
                p_background *= 0.25

            # Corrected probability for background probability
            probabilities.append(p_motif / p_background)

        # Choose a random position based on probability distribution
        sequence['motif_position'] = self.choose_random(probabilities)

        return


    def calculate_entropy(self, motif):

        """
        Calculate the relative entropy of the given motif matrix.
        """

        entropy = 0

        for base in "ATCG":
            for i in range(self.motif_width):
                entropy -= motif[base][i] * math.log(motif[base][i], 2)

                # I think there'se no need to use the background pseudocounts
                # here, because we already used them in the motif matrix.
                # entropy += motif[base][i] * log(motif[base][i] /
                #                                 self.pseudocounts[base], 2)

        return entropy


    def choose_random(self, distribution):

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
