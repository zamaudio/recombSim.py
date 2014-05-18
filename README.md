recombsim
=========

This program builds a random DNA sequence as an ancestral sequence.
It then generates 3 recombinant subsequences and replaces the ancestral sequence
at known locations.  Then by a random process, this program will generate progressively
noisier sequences by mutating at random locations, from 1% up to 100% of the sequence.

It outputs the ancestral sequence in one file, and an alignment of the error prone
sequences in order of randomness in a second file.

The resulting dataset can be used to test the accuracy of various algorithms for
detecting homologous recombination events by checking if the results match the
original regions of inserted recombinant subsequences.
Also, the sensitivity can be tested by observing at which point the algorithm
breaks down and no longer is able to predict the recombination events.

Dependencies:
=============
Bio++ (libbpp-core, libbpp-seq) found here: http://biopp.univ-montp2.fr/

A C++ compiler


The pseudo random number generator used to randomise the sequence is Mersenne Twister.
