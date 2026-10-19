# BINF6250 - Project 05

# Introduction
This project implements the Smith-Waterman algorithm, a dynamic programming method for local pairwise sequence alignment that finds the highest-scoring shared region between two sequences under a customizable scoring system (defaults: match +1, mismatch -1, gap -1). The program fills a scoring matrix cell by cell, where each cell takes the largest of four options (up-left neighbor plus match or mismatch, up neighbor plus gap, left neighbor plus gap, or 0, which keeps scores from going negative), records the winning direction in a second matrix, and rebuilds the alignment by tracing those directions back from the highest score until it reaches a 0. The work is split across three functions: `cal_score()` scores one cell and returns its direction code (0 = END, 1 = DIAG, 2 = UP, 3 = LEFT), `traceback()` follows the directions to build the two aligned sequences, and `smith_waterman()` builds both NumPy matrices, fills them, finds and prints the best score, calls `traceback()`, and returns the results.  All required functionality described in the assignment instructions was completed. All three functions were implemented, and each was tested with a sanity check on a small example worked out by hand before the full program was run. Test 1 matches the expected output provided in the notebook, and Test 2, which has no expected output, was checked by reasoning through the scores in the printed matrix.

# Pseudocode
### cal_score()
```
Get this row's letter: seq1 at position (i - 1)
Get this column's letter: seq2 at position (j - 1)

Diagonal number:
    IF this row's letter equals this column's letter:
        diag_score = matrix at (i - 1, j - 1) + match
    ELSE:
        diag_score = matrix at (i - 1, j - 1) + mismatch
up_score = matrix at (i - 1, j) + gap
left_score = matrix at (i, j - 1) + gap

score = biggest of: 0, diag_score, up_score, left_score

Pick the move (check in this order, first match wins if tied):
    IF score is 0:          move = 0 (END)
    ELIF score = diag:      move = 1 (DIAG)
    ELIF score = up:        move = 2 (UP)
    ELSE:                   move = 3 (LEFT)

Return score and move
```

### traceback()
```
aligned_seq1 = empty string
aligned_seq2 = empty string
current_row, current_col = maximum_position
current_move = traceback_matrix at (current_row, current_col)

WHILE current_move is not 0 (END):
    IF current_move is 1 (DIAG):
        add seq1 letter at (current_row - 1) to front of aligned_seq1
        add seq2 letter at (current_col - 1) to front of aligned_seq2
        move up one row and left one column
    ELIF current_move is 2 (UP):
        add seq1 letter at (current_row - 1) to front of aligned_seq1
        add "-" to front of aligned_seq2
        move up one row
    ELIF current_move is 3 (LEFT):
        add "-" to front of aligned_seq1
        add seq2 letter at (current_col - 1) to front of aligned_seq2
        move left one column

    current_move = traceback_matrix at (current_row, current_col)

Return aligned_seq1, aligned_seq2
```

### smith_waterman()
```
rows = length of seq1 + 1
cols = length of seq2 + 1
score_matrix = grid of zeros, size rows x cols (holds decimals)
traceback_matrix = grid of zeros, size rows x cols (whole numbers 0, 1, 2, 3)

FOR i from 1 to rows - 1:
    FOR j from 1 to cols - 1:
        score, move = cal_score(score_matrix, seq1, seq2, i, j, match, mismatch, gap)
        score_matrix at (i, j) = score
        traceback_matrix at (i, j) = move

best_score = highest number in score_matrix
max_pos = (row, column) where best_score is (if tied, keeps only the first one found)

IF best_score is 0: print that no local alignment was found

print best_score

aligned_seq1, aligned_seq2 = traceback(seq1, seq2, traceback_matrix, max_pos)

Return aligned_seq1, aligned_seq2, score_matrix
```


# Results Output
### Test 1: TACTTAG vs CACATTAA (Example from slides)
```
Alignment score:  4.0
AC-TTA
ACATTA
[[0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 1. 1. 0. 0.]
 [0. 0. 1. 0. 1. 0. 0. 2. 1.]
 [0. 1. 0. 2. 1. 0. 0. 1. 1.]
 [0. 0. 0. 1. 1. 2. 1. 0. 0.]
 [0. 0. 0. 0. 0. 2. 3. 2. 1.]
 [0. 0. 1. 0. 1. 1. 2. 4. 3.]
 [0. 0. 0. 0. 0. 0. 1. 3. 3.]]
```

The alignment and every matrix value match the expected output provided in the notebook. The matrix prints as `0.` rather than `0` because the scoring matrix holds decimals, but the values are identical.

### Test 2: GGGGAAAAAAAAAAAACCCC vs GGGGAAAAAACCCC
```
Alignment score:  10.0
GGGGAAAAAA
GGGGAAAAAA
[[ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
 [ 0.  1.  1.  1.  1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
 [ 0.  1.  2.  2.  2.  1.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
 [ 0.  1.  2.  3.  3.  2.  1.  0.  0.  0.  0.  0.  0.  0.  0.]
 [ 0.  1.  2.  3.  4.  3.  2.  1.  0.  0.  0.  0.  0.  0.  0.]
 [ 0.  0.  1.  2.  3.  5.  4.  3.  2.  1.  1.  0.  0.  0.  0.]
 [ 0.  0.  0.  1.  2.  4.  6.  5.  4.  3.  2.  1.  0.  0.  0.]
 [ 0.  0.  0.  0.  1.  3.  5.  7.  6.  5.  4.  3.  2.  1.  0.]
 [ 0.  0.  0.  0.  0.  2.  4.  6.  8.  7.  6.  5.  4.  3.  2.]
 [ 0.  0.  0.  0.  0.  1.  3.  5.  7.  9.  8.  7.  6.  5.  4.]
 [ 0.  0.  0.  0.  0.  1.  2.  4.  6.  8. 10.  9.  8.  7.  6.]
 [ 0.  0.  0.  0.  0.  1.  2.  3.  5.  7.  9.  9.  8.  7.  6.]
 [ 0.  0.  0.  0.  0.  1.  2.  3.  4.  6.  8.  8.  8.  7.  6.]
 [ 0.  0.  0.  0.  0.  1.  2.  3.  4.  5.  7.  7.  7.  7.  6.]
 [ 0.  0.  0.  0.  0.  1.  2.  3.  4.  5.  6.  6.  6.  6.  6.]
 [ 0.  0.  0.  0.  0.  1.  2.  3.  4.  5.  6.  5.  5.  5.  5.]
 [ 0.  0.  0.  0.  0.  1.  2.  3.  4.  5.  6.  5.  4.  4.  4.]
 [ 0.  0.  0.  0.  0.  0.  1.  2.  3.  4.  5.  7.  6.  5.  5.]
 [ 0.  0.  0.  0.  0.  0.  0.  1.  2.  3.  4.  6.  8.  7.  6.]
 [ 0.  0.  0.  0.  0.  0.  0.  0.  1.  2.  3.  5.  7.  9.  8.]
 [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  1.  2.  4.  6.  8. 10.]]
```

The printed matrix contains two cells with the best score of 10: one at row 10 / column 10, which traces back to `GGGGAAAAAA`, and one in the bottom-right corner, which would trace back to `AAAAAACCCC`.

# Successes

# Struggles

# Reflection and Adaptive Learning

# Personal Reflections
## Group leader

## Other member

## Other member: Ildiko Polyak
Compared to graph algorithms, Smith-Waterman felt much easier to understand and code, since each cell of the matrix is just three additions and a max. My strategy has been to put most of my effort into understanding the algorithm first and then writing robust pseudocode, so that writing the real code (my biggest weakness) is easier translation. Where I struggled this week was not the algorithm itself but how the project split it across three functions, because the steps were not divided in the same order I had learned them. What helped was writing out a generic pseudocode for Smith-Waterman pairwise alignment first, in its natural order (set up, fill, find the best cell, trace back, output), and then assigning each step to the function it belonged in, which kept clear which steps needed to be written where and in what order. For the same reason, I added a simple print test after each function, running it on a small example I had already worked out by hand, so I could confirm each function worked on its own instead of waiting to troubleshoot the entire program at the end. Reflecting on that, I can still only process an algorithm in its proper order, and trying to separate it into pieces for a more "advanced" program leads to confusion for me. In the beginning of the semester, this stressed me out because I felt like I had to "master" each algorithm, but the pace of the class is fast and switches to a new world the next week.  Halfway throughout the semester, I'm now OK with understanding the algorithm at my level and continued use and consideration will fortify it.  I also need to spend more time learning how to work with numpy.

Even though the group's schedule had very little overlap for much collaborative work, we were able to meet for two short sessions and it was nice to reflect on our experience in the class together and not feel so alone in the fast pace and struggle to keep up with new difficult information every week.  It helped build confidence in the learning process without the need to master an outcome.

# Generative AI Appendix
Anthropic Claude Opus 5.5 was used for this project.

Ildiko is inexperienced with NumPy and used Claude to write the NumPy-specific code, since she did not yet know those functions: creating the scoring and traceback matrices (`np.zeros()`, including `dtype=int` for the traceback matrix), finding the best score and its position (`np.max()`, `np.argmax()`, and `np.unravel_index()`), and building the small test matrix for the `traceback()` sanity check (`np.array()`). Claude also explained how each of these functions works before it was added to the code.
