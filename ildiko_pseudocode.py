#General Smith-Waterman psuedocode framework for project pseudocode

#Inputs:
    #sequence A (goes down the left side)
    #sequence B (goes across the top)
    #match score, mismatch score, gap penalty

#Step 1: Set up - smith_waterman()
    #Make a score grid with (length of A + 1) rows and (length of B + 1) columns
    #Fill every cell with 0
    #Make an second grid for arrows the same size, fill every cell with "END"

#Step 2: Fill 
    #FOR each row, from row 1 to last row:
        #FOR each column, from column 1 to last column: - smith_waterman()
            #Above-left number = above-left neighbor's score + (match score if this row's letter equals this column's letter, otherwise mismatch score) - cal_score()
            #Above number = above neighbor's score + gap penalty
            #Left number = left neighbor's score + gap penalty
            #Zero = 0 (resets so there is no negative on the grid)
            #This cell's score = biggest number from the 4 options above

            #Record this cell's arrow: 
                #If this cell's score is 0             END
                #Otherwise, if above-left number won   ↘ 
                #Otherwise, if above number won        ↓ 
                #Otherwise, if left number won         →

#Step 3: Find the best cell - smith_waterman()
    #Look at every cell in the score grid
    #find the highest score
    #Remember that score and its row and column (this is where traceback starts)

#Step 4: Trace back - traceback()
    #Start two empty aligned strings: aligned A and aligned B
    #Start at the best cell's row and column
    #WHILE the current cell's arrouw is not END:
        #IF the arrow is up-left:
            #Add this row's letter (from A) to the front of aligned A
            #Add this column's letter (from B) to the front of aligned B
            #Move up one row and left one column
        #IF the arrow is up:
            #Add this row's letter (from A) to the front of aligned A
            #Add a gap "-" to the front of aligned B
            #Move up one row
        #If the arrow is left:
            #Add a gap "-" to the front of aligned A
            #Add this column's letter (from B) to the front of aligned B
            #Move left one column

#Step 5:  Output - smith_waterman()
    #Return aligned A, aligned B, the best score, and the score grid


#Pseudocode specific for Project 05
import numpy as np

def cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap): #Step 2: Fill, scoring one cell (four numbers, biggest, arrow)
    '''Calculate score for position (i,j) in scoring matrix, also record move to trace back
    
    Args:
        matrix (numpy array): scoring matrix
        seq1 (str): sequence 1
        seq2 (str): sequence 2
        i (int): current row number
        j (int): current column number
        
    Returns:
        score in position (i,j)    
        move to trace back: 0-END, 1-DIAG, 2-UP, 3-LEFT
        
    Pseudocode:
        Calculate scores based on upper-left, up, and left neighbors:
            diag_score = upper-left + (match or mismatch)
            up_score = up + gap
            left_score = left + gap
        score = max(0, diag_score, up_score, left_score)
        traceback = maximum direction or end
        
    '''
    #Get this row's letter: seq1 at position (i-1)
    #Get this column's letter: seq2 at position (j-1)
    
    #Diagonal number:
        #If this row's letter equals this column's letter:
            #diag_score = matrix at (i-1, j-1) + match
        #Otherwise:
            #diag_score = matrix at (i-1, j-1) + mismatch
    #up_score = matrix at (i-1, j) + gap
    #left_score = matrix at (i, j-1) + gap

    #score = biggest of: 0, diag_score, up_score, left_score

    #Pick the move (check in this order, first match wins if tie score):
        #IF score is 0:  move = 0 (END)
        #ELIF score = diag:  move = 1 (DIAG)
        #ElIF score = up:  move = 2 (UP)
        #ELSE:  move = 3 (LEFT)

    #Return score and move

def traceback(seq1, seq2, traceback_matrix, maximum_position): #Step 4: Trace back
    '''Find the optimal path through scoring marix
        
        diagonal: match/mismatch
        up: gap in seq1 #LANDMINE - GAP IN SEQ2
        left: gap in seq2 #LANDMINE - GAP IN SEQ1
        
    Args:
        seq1 (str) : First sequence being aligned
        seq2 (str) : Second sequence being aligned
        traceback_matrix (numpy array): traceback matrix
        maximum_position (tuple): starting position to trace back from
        
    Returns:
        aligned_seq1 (str): e.g. GTTGAC
        aligned_seq2 (str): e.g. GTT-AC
        
    Pseudocode:
        while current_move != END:
            current_move = traceback_matrix[current_row][current_col]
            if current_move == DIAG:
                ...
            elif current_move == UP:
                ...
            elif current_move == LEFT:
                ...
            
    '''
    #aligned_seq1 = empty string
    #aligned_seq2 = empty string
    #current_row, current_col = maximum_position (row, column)
    #current_move = traceback_matrix at (current_row, current_col)

    #WHILE current_move is not 0 (END):
        #IF current_move is 1 (DIAG):
            #add seq1 letter at (current_row - 1) to front of aligned_seq1
            #add seq2 letter at (current_col - 1) to front of aligned_seq2
            #move up one row and left one column

        #ELIF current_move is 2 (UP):
            #add seq1 letter at (currnet_row - 1) to front of aligned_seq1
            #add "-" to front of aligned_seq2
            #move up one row

        #ELIF current_move is 3 (LEFT): 
            #add "-" to front of aligned_seq1
            #add seq2 letter at (current_col - 1) to front of aligned_seq2
            #move left one column

        #current_move = traceback_matrix at (current_row, current_col)

    #Return aligned_seq1, aligned_seq2


def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-1): #Step 1: Set up grids; Step 2: Fill, the loop through every cell; Step 3: Find best cell; Step 5: Output
    '''Smith-Waterman algorithm for local alignment
    
    Args:
        seq1 (str): input seq 1
        seq2 (str): input seq 2
        match: default = +1
        mismatch: default = -1
        gap: default = -1
    
    Returns:
        aligned_seq1 (str)
        aligned_seq2 (str)
        score_matrix (numpy array): scoring matrix
    '''
    #rows = length of seq1 + 1
    #cols = length of seq2 + 1
    #score_matrix = grid of zeros, size rows x cols (decimals, not whole number)
    #traceback_matrix = grid of zeros, size rows x cols (whole numbers 0, 1, 2, 3)

    #for i from 1 to rows -1:
        #for j from 1 to cols -1:
            #score, move = cal_score (score_matrix, seq1, seq2, i, j, match, mismatch, gap)
            #score_matrix at (i, j) = score
            #traceback_matrix at (i, j) = move

    #best_score = highest number in score matrix
    #if best_score is 0 (END) display that and makes sure it doesn't get caught up in loop
    #max_pos = (row, column) where best_score is (if there is a tie, this keeps only the first one found)

    #print best_score

    #aligned_seq1, aligned_seq2, = traceback (seq1, seq2, traceback_matrix, max_pos)

    #Return aligned_seq1, aligned_seq2, score_matrix

# Example from slides
seq1 = 'TACTTAG'
seq2 = 'CACATTAA'

aligned_seq1, aligned_seq2, score_matrix = smith_waterman(seq1,seq2)

print (aligned_seq1)
print (aligned_seq2)
print (score_matrix)

seq3 = "GGGGAAAAAAAAAAAACCCC"
seq4 = "GGGGAAAAAACCCC"

aligned_seq3, aligned_seq4, score_matrix = smith_waterman(seq3,seq4)

print (aligned_seq3)
print (aligned_seq4)
print (score_matrix)
