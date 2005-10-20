"""
Set of functions implementing the Gibbs sampling algorithm to find local
alignments of sequences.

TODO: more text here

Martijn Vermaat mvermaat@cs.vu.nl
"""


import random
import math


PS_ITERATIONS_FACTOR = 4
ONE_PLUS_ONE = 2


class Gibbs:


    """
    This class implements the Gibbs sampling algorithm.
    """


    __sequences = []
    __motif_width = 0
    __pseudocounts = {}


    def __init__(self, sequences, motif_width, pseudocounts_weight=0.1):

        """
        Gibbs takes a list of sequence objects, the width of the motif to find
        and the weight to use for pseudocounts (default 0.1).

        The sequence object are dictionaries with keys:

            title           title of the sequence (e.g. gene name)
            sequence        string with ATCG occurences in sequence
            motif_position  integer position of the motif in sequence
        """

        self.__sequences = sequences
        self.__motif_width = motif_width

        self.__calculate_pseudocounts(pseudocounts_weight)

        return


    def find_motif(self, iterations=50, phase_shifts=0, ps_frequency=12):

        """
        Finds a motif in sequences using Gibbs sampling.

        Parameters:

            iterations    non-improving iterations before quiting (50)
            phase_shifts  number of phase shifts to examine (0)
            ps_frequency  frequency of iterations detecting phase shifts (12)

        Stores the best alignment found.
        """

        #self.__random_motif_positions()

        self.__initial_motif_positions()

        best_entropy = 0
        best_alignment = [0] * len(self.__sequences)

        # TODO: this is only for graphing
        entropies = []

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

            #print_motif(motif)

            entropies.append(self.__calculate_entropy(motif))

            self.__calculate_position(motif, current_sequence)

            # Do a phase shift every ps_frequency iterations
            if phase_shifts > 0 and j % ps_frequency == 0:
                shift = self.__phase_shift(phase_shifts)
                # Give the algorithm a chance of recovering from
                # badly choosen phase shifts
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

        motif = self.__calculate_motif(self.__sequences)
        entropy = self.__calculate_entropy(motif)

        #self.__print_entropies(entropies)

        print best_entropy

        self.__print_motif(motif)

        self.__print_sequences()

        return


    def __phase_shift(self, shift):

        """
        Calculate motif for some small phase shifts in each sequence and apply
        the fase shift with best entropy.
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

        #print "ja: %i (%i, %i)" % (best - shift_left, shift_left, shift_right)

        return best - shift_left


    def __calculate_pseudocounts(self, weight):

        """
        Calculate for each base a weighted pseudocount.
        """

        # Total number of bases in sequences
        total = float(sum([len(s['sequence']) for s in self.__sequences]))

        for base in "ATCG":
            # Number of occurences of this base in sequences
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


    def __initial_motif_positions(self):

        """
        Populate the list of sequences with a random position of the motif for
        each sequence based on the following heuristic.

        The idea is that significant motifs are mostly the ones that contain
        relatively a lot of bases with a low background frequency.
        In order to guide the algorithm to regions in the sequences where more
        of these bases occur, we do the following:

            Pick the base with the lowest background frequency. Now, in each
            sequence, find all words of length pattern_width with at least
            number_of_occurences occurences of that base. Choose a random one
            of these position and initialize the motif position of the
            sequence by aligning it with the random choosen position.

        Of course, significant motifs might not always contain a certain
        number of occurences of the base with lowest background frequency, let
        alone within a certain amount of bases apart from each other.
        """

        # TODO: command line options for these
        pattern_width = 6
        number_of_occurrences = 2

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

                    # Count occurences of lowest_base in word
                    n = s['sequence'][r : r+pattern_width].count(lowest_base)

                    if n >= occurrences:
                        positions.append(r)

                if len(positions) > 0:
                    break

                # Try the same with one occurence of lowest_base less
                occurrences -= 1

            if len(positions) > 0:

                position = random.choice(positions)
                # TODO: center motif and word here, look out for end of
                # sequence

            else:

                # This will not happen a lot, but choose a random position if
                # the sequence doesn't have at least one occurence of
                # lowest_base
                position = random.randint(
                    0, len(s['sequence']) - self.__motif_width)

            s['motif_position'] = position

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

        # TODO: remove this or make it a debug thing

        from rpy import r

        r.postscript("entropies.ps")
        r.plot(entropies, type='b', xlab="Iterations", ylab="Entropy")
        r.dev_off()

        return


    def __print_motif(self, motif):

        # Print position weight matrix for motif
        for base in "ATCG":
            print base,
            # Do not uncomment this comment
            for weight in motif[base]:
                print "%1.2f" % weight,   # we have floats
            print

        return


    def __print_sequences(self):

        # Print motif occurence in each sequence and motif position
        for s in self.__sequences:
            start, end = s['motif_position'], s['motif_position'] + self.__motif_width
            print "%s... motif at position %s" % (s['sequence'][start:end], s['motif_position'])


class GibbsError(Exception):
    """
    Exception thrown by Gibbs instances.
    """
    # Our own exception doesn't need anything fancy.
    pass
