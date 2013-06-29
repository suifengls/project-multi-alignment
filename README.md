Project: multi-alignment
=====

We discussed a space-saving technique for pairwise sequence alignment
based on divide-and-conquer last week. The technique can be easily
extended to multiple sequence alignment. Your task is to design an
O(m\*n) space, O(m\*n\*k) time SP alignment algorithm for three DNA
sequences of lengths m, n, and k, respectively. 

(i)  For simplicity, let's consider global alignment and the basic 
SP score model where gaps are not specially treated. 

(ii) Your program should work for any score function on nucleotides.
In other words, the user should be able to input a score function
in the form of a 5x5 matrix indexed by A, C, G, T, and space.
The SP-score of a column of letters/spaces is the sum of the scores
of each pair of letters/spaces in the column.

To test your program, use the Blast scores: Match = 5, Mismatch = -4,
and Indel = -8. The score between two spaces is 0.

(iii) Test your program on the following five sets of sequences 
(500 bps - 6000 pbs)
