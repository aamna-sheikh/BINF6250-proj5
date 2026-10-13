## A description of the project


## Pertinent documentation and usage write-up


## Pseudocode


## Reflection narrative that describes what went right, what went wrong, what didn't you understand, how it could have been better, and how you could have been better.
### What went right

### What went wrong
The recursive walk worked on the small examples but crashed on the yeast data with RecursionError, because Python only allows about 1,000 nested calls and the yeast walk needed millions.
I fixed it by rewriting the walk as a loop with a stack (a list I add to and take from). My teammate fixed it differently, by raising the recursion limit to 25 million.
The yeast data has 30x coverage, so every k-mer appears about 30 times. Without checking for duplicates in add_edge, the graph would have 30 copies of every edge.

My teammate and I got very different results (30 contigs vs 7,632) from the same data, and at first we didn’t know why or which one was right.
We found that our graphs were identical (the total lengths differ by exactly 30 letters per extra contig), and the only difference was one word: if vs while in the walk.
My result has one contig of 22.6 million letters, which is longer than the whole genome. This showed me that a correct Eulerian walk is not always a realistic genome.

### Things I still don't understand
How real assemblers decide where to stop a contig at a fork, and how they remove the reverse-strand duplicates.

## Appendix (AI usage, citations, etc.)

