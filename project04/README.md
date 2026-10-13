## A description of the project
## Description

In this project, we built a genome assembler using a **de Bruijn graph**.
A sequencing machine reads a genome in millions of short pieces called **reads**, in random order. Assembly means putting those pieces back together to rebuild the original DNA.

Our program does this in four steps:
1. **Read the data:** load the DNA reads from a FASTQ file.
2. **Build the graph:** cut every read into small overlapping pieces of length k (called k-mers). Each k-mer is split into a left half and a right half, and we draw an arrow from the left half to the right half. Reads that overlap share the same pieces, so they connect automatically.
3. **Walk the graph:** follow the arrows, using each arrow exactly once (an Eulerian walk).
4. **Rebuild the DNA:** turn each walk back into a DNA sequence, called a **contig**.

We first tested the program on small examples, then ran it on 2.4 million perfect 150 bp reads from the yeast (*S. cerevisiae*) genome, using k = 31. The program writes the contigs to a FASTA file and saves a statistics report (number of contigs, total length, N50, and more).
We assumed perfect sequencing, which means no errors and no variants in the reads.

## Pertinent documentation and usage write-up
### Requirements
- Python 3
- Jupyter Notebook (or an editor that runs notebook cells, like VS Code or PyCharm)
- No extra packages are needed. The code only uses Python's built-in modules: `gzip`, `collections`, `random`, `time`, and `datetime`.

### Files
- `project04_yeast.ipynb`: the notebook with all the code
- `data/sacCer_R64_150SE_30x_wgsim.fq.gz`: the yeast reads (from `project04.zip`)

### How to run
1. Unzip `project04.zip` into the project folder, so the reads are in the `data/` folder.
2. Open the notebook.
3. Run the cells **in order, from top to bottom**:
   - the `read_fastq` cell
   - the `DeBruijnGraph` class cell
   - the toy example cell (a quick check that should give 1 contig)
   - the real data (yeast) driver cell

### Settings
These are at the top of the driver cell:
Setting- Value we used (Meaning)
`INPUT_FILE`- `data/sacCer_R64_150SE_30x_wgsim.fq.gz` (the reads file)
`K_MER_SIZE`- `31` (the length of each k-mer)
`GENOME_SIZE`- `12,071,326` (the size of the yeast genome, used to compare results)
`OUTPUT_FASTA`- `genome_assembly.fasta` (where the contigs are saved)
`STATS_FILE`- `genome_assembly_stats.txt` (where the statistics are saved)

### Output
- `genome_assembly.fasta`: every contig, each with a header line like `>contig_1 length=22639116`
- `genome_assembly_stats.txt`: a report with the number of contigs, total length, N50, longest contig, and timing

### Run time
On our computer, the yeast run took about 5 minutes:
- loading reads: about 2 seconds
- building the graph: about 3 minutes
- walking the graph: about 2 minutes

It also uses a lot of memory, because the graph has about 23 million nodes.

### What each method does

Method -What it does 

`read_fastq` -reads the FASTQ file and keeps only the DNA lines |
`add_edge` -adds an arrow between two nodes (only once, no duplicates) |
`remove_edge` -removes an arrow |
`build_graph_from_reads` -cuts every read into k-mers and adds an arrow for each one |
`eulerian_walk` -walks the graph using each arrow once, using a loop and a stack |
`assemble_contigs` -finds starting points, walks from each one, and turns each walk into a contig |
`tour_to_sequence` -turns a list of nodes back into a DNA string |
`get_assembly_stats` -calculates the number of contigs, lengths, mean, and N50 |
`write_fasta` -saves the contigs to a FASTA file |

### Small example
```python
reads = ["ATGGC", "GGCGT"]
dbg = DeBruijnGraph(reads, k=4)
contigs = dbg.assemble_contigs(seed=42)
print(contigs)   # ['ATGGCGT']
```

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

