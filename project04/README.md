## A description of the project
## Description

In this project, we built a genome assembler using a **de Bruijn graph**.
A sequencing machine can't read a whole genome from start to finish. Instead, it reads millions of short pieces, called **reads**, from random places in the genome. Assembly means putting those pieces back together to rebuild the original DNA.

Our program does this in six steps:
1. **Read the data:** load the DNA reads from a FASTQ file.
2. **Build the graph:** cut every read into small overlapping pieces of length *k*, called **k-mers**. Each k-mer is split into a left part (its first *k*−1 letters) and a right part (its last *k*−1 letters), and we draw an edge from the left part to the right part. The parts are the **nodes**. Reads that overlap share the same pieces, so they connect automatically.
3. **Walk the graph:** follow the edges, using each edge exactly once (an Eulerian walk).
4. **Rebuild the DNA:** turn each walk back into a DNA sequence, called a **contig**. Steps 3 and 4 repeat until every edge in the graph has been used, so every part of the graph ends up in a contig.
5. Remove reverse-strand copies: DNA has two strands, and our reads come from both, so most regions get assembled twice. We keep one copy of each and drop the other.
6. Measure and save: calculate statistics about the contigs (number of contigs, lengths, N50) and save the contigs to a FASTA file.
   
We first tested the program on small examples, then ran it on 2,414,266 simulated 150 bp reads from the baker's yeast (*Saccharomyces cerevisiae*) R64 genome, which is 12,071,326 bp long, using k = 31. That is about 30× coverage, meaning each position was read about 30 times on average. The program writes the contigs to a FASTA file and saves a statistics report (number of contigs, total length, N50, and more).

We assumed perfect sequencing: no sequencing errors and no variants in the reads.

## Pertinent documentation and usage write-up
### Requirements
- Python 3.11 or newer (we ran it on Python 3.14)
- Jupyter Notebook (or an editor that runs notebook cells, like VS Code or PyCharm)
- No extra packages are needed. The code only uses Python's built-in modules: `gzip`, `collections`, `random`, `time`, and `datetime`.

### Usage

Open `project04_yeast.ipynb` in Jupyter and run all cells in order. The notebook is split into:

- Optional check: **Peek at the data** - prints the first 2 reads of the FASTQ file, to confirm the file is in the right place and in the expected format.\
- Core deliverables: **`read_fastq()`** - reads the FASTQ file, and the **`DeBruijnGraph`** class - builds the graph, walks it, assembles the contigs, removes reverse-strand copies, and calculates statistics.\
- Test: **Toy example** - runs `DeBruijnGraph` on 9 small overlapping reads. It should print one contig, `ATGGCGTACGTTACCATGGGCCTAAA` (26 bp), and `Correct: True`.\
- Provided program: **Driver Program** - runs `DeBruijnGraph` on the yeast reads, prints the results, and saves the output files. We changed one line: after `contigs = dbg.assemble_contigs()`, we added `contigs = dbg.remove_reverse_complement_duplicates(contigs)`.\
- Optional check: **Peak memory** - prints the most memory the notebook used during the run.\

Input data files are NOT tracked in this repo. Unzip `project04.zip` and place the reads in a `data/` folder next to the notebook before running it, as shown in Project structure below.\

### Project structure

```
project04/
├── project04_yeast.ipynb
├── README.md
└── data/                                   (not tracked; from project04.zip)
    └── sacCer_R64_150SE_30x_wgsim.fq.gz
```

Running the driver creates `genome_assembly.fasta` and `genome_assembly_stats.txt` in the same folder as the notebook.

### Settings
These are at the top of the driver cell in the notenook:
Setting- Value we used (Meaning)
`INPUT_FILE`- `data/sacCer_R64_150SE_30x_wgsim.fq.gz` (the reads file)
`K_MER_SIZE`- `31` (the length of each k-mer)
`GENOME_SIZE`- `12,071,326` (the size of the yeast genome, used to compare results)
`OUTPUT_FASTA`- `genome_assembly.fasta` (where the contigs are saved)
`STATS_FILE`- `genome_assembly_stats.txt` (where the statistics are saved)

### Output
- `genome_assembly.fasta`: every contig, longest first, each with a header line like `>contig_1 length=233925`, followed by its sequence
- `genome_assembly_stats.txt`: a report with the number of contigs, total length, N50, longest contig, and timing

### Run time and memory

Times from the final run of `project04_yeast.ipynb` on a 16 GB Mac:

| Step | Time |
|---|---|
| Load reads | 3.84 s |
| Build graph | 313.09 s |
| Assemble (including removing reverse-complement duplicates) | 96.99 s |
| **Total** | **413.92 s (6.90 min)** |
| **Peak memory** | **4.74 GB** (measured with `resource.getrusage`) |

The graph has 22,950,037 nodes and 22,957,882 edges, which is why it needs several GB of memory. Your times will depend on your computer.

### What each method does

As written in `project04_yeast.ipynb`:

| Method | What it does |
|---|---|
| `read_fastq` | Reads the FASTQ file and keeps only the DNA line of each read (line 2 of every 4). Works with plain or `.gz` files. |
| `__init__` | Makes an empty graph, saves k, and calls `build_graph_from_reads`. |
| `add_edge` | Adds an edge from one node to another. |
| `remove_edge` | Removes one copy of an edge, if it exists. Implemented and tested, but not called by the walk, which removes edges from the graph it is given directly. |
| `build_graph_from_reads` | Cuts every read into k-mers and adds an edge for each one, skipping edges that are already in the graph. |
| `eulerian_walk` | Recursive. If the current node has an unused edge, picks one at random, removes it, and calls itself on the next node. Then adds the current node to the tour. Every call writes into one shared list (`backward_tour`), so nothing is copied. Returns the nodes in reverse order. |
| `tour_to_sequence` | Turns a list of nodes back into a DNA string: the first node whole, then the last letter of each node after it. |
| `assemble_contigs` | Counts edges coming into each node, raises the recursion limit to fit the graph (and puts it back afterward), then walks from nodes with more edges out than in, and then from any node with edges left, until every edge is used. Flips each walk and turns it into a contig. By default it walks `self.graph` directly to save memory, which empties it; `copy_graph=True` walks a copy instead. |
| `remove_reverse_complement_duplicates` | Goes through the contigs longest first and drops any contig where at least 90% of its k-mers are reverse complements of k-mers already kept, so each region is kept once instead of once per strand. |
| `get_assembly_stats` | Calculates the number of contigs, total, longest, shortest and mean length, and N50. |
| `write_fasta` | Saves the contigs to a FASTA file. |

### Small example
```python
reads = ["ATGGC", "GGCGT"]
dbg = DeBruijnGraph(reads, k=4)
contigs = dbg.assemble_contigs(seed=42)
print(contigs)   # ['ATGGCGT']
```
The k-mers are `ATGG`, `TGGC`, `GGCG` and `GCGT`, which give the edges `ATG → TGG → GGC → GCG → CGT`. Walking that path and spelling it out gives `ATGGCGT`.

## Pseudocode
Check the data file (optional)
```
    Open the compressed FASTQ file as text
    Print the first 8 lines, which are the first 2 reads (4 lines each):
        line 1: the read's name, starting with @
        line 2: the DNA sequence
        line 3: a + sign
        line 4: the quality scores
```
Read the data
```
read_fastq(file):
    for every 4 lines in the file:
        keep line 2 (the DNA), uppercase, without the newline
    return the list of reads
```
De Bruijn Graph Genome Assembly Module:
```
# Set up map and Helpers functions
__init__( reads, k):
     Make an empty graph
        Save k
        Build the graph from the reads (below)
add_edge(left, right): 
     Add an edge from left to right
        Put right on left's list of next nodes
remove_edge( left, right): #Written and tested, but not called by the walk, which removes edges from the graph it is given directly
     Remove an edge from left to right               
        If right is on left's list:
             Take one copy of it off

# Build the graph
build_graph_from_reads(reads, k):
    for each read:
        for each k-mer in the read (k letters, sliding over by one letter each time):
            left node  = first k-1 letters (the k-mer without its last letter)
            right node = last k-1 letters (the k-mer without its first letter)
            If there is no edge from left to right yet: Add one 
                add_edge(left, right)

# Walk the graph
eulerian_walk(node, graph, backward_tour):
    if backward_tour is not given:
        backward_tour = empty list # first call makes the shared list
    if node has an unused edge:
        next_node = pick one at random
        remove that edge
        eulerian_walk(next_node, graph, backward_tour) # recurse
    add node to backward_tour
    return backward_tour  # in reverse order

# Rebuild the DNA
tour_to_sequence(tour):
    sequence = first node
    for each node after the first:
        add its last letter
    return sequence

# Assemble all the contigs
assemble_contigs(copy_graph = False):
    graph = a copy of the graph if copy_graph, otherwise the graph itself
    count edges coming into each node
    raise the recursion limit to (number of edges + 1000)
    try:
        start nodes = nodes with more edges out than in, then every node
        for each start node:
            while it still has an unused edge:
                tour = eulerian_walk(start), flipped to forward order
                add tour_to_sequence(tour) to the contigs
    finally:
        put the recursion limit back
    return the contigs

# Remove reverse-strand copies
remove_reverse_complement_duplicates(contigs, cut off = 0.9):
    seen = empty set of k-mers
    for each contig, longest first:
        flip it to the other strand (swap A/T and C/G, then reverse)
        if at least 90% of the flipped contig's k-mers are in seen:
            skip it                           # it's the other strand's copy
        else:
            keep it, and add its k-mers to seen
    return the kept contigs

# Measure and save the result
get_assembly_stats(contigs):
    sort the contig lengths from longest to shortest
    N50 = the length at which the running total reaches half the total length
    return the counts, lengths, mean, and N50
write_fasta(contigs, file):
    for each contig, numbered from 1:
        write ">contig_N length=L", then the sequence
```

## Reflection narrative that describes what went right, what went wrong, what didn't you understand, how it could have been better, and how you could have been better.
### What went right
Testing every function on small, hand-checkable examples before running the yeast data caught most of our bugs early, when they were easy to understand. The toy example (9 reads that should give one 26 bp contig) became our 1st check: if it didn't give exactly that contig, something was wrong.
Comparing our three notebooks was also useful. Because we each got different results from the same data, we had the opportunity to dig and work out why, which is helpful. 

### What went wrong
**Duplicate edges.** The yeast data has about 30× coverage, so most k-mers appear about 30 times. Without a duplicate check, the graph would hold about 30 copies of most edges, and the walk would drive over the same stretches again and again. On the toy example, this gave a 50 bp nonsense contig instead of the correct 26 bp one.

**The recursion limit.** The recursive walk worked on the small examples but crashed on the yeast data with `RecursionError`, because Python only allows about 1,000 nested calls by default. We fixed this in different ways: 
- One of us rewrote the walk as a loop with a stack (a list we add to and take from),
- Another kept the recursion but raised the limit to 25mil,
- The third also kept the recursion, but set the limit from the number of edges in the graph. Each recursive call uses one edge, so the walk can never need more calls than there are edges. This version also puts the old limit back afterward. We chose this approach for our final notebook, and added try/finally so the limit is restored even if the walk crashes. It would fit any data set, from the toy example to the full yeast genome, without guessing a number, and it doesn’t leave Python’s setting changed after the program finishes.

**Memory.** The graph for the yeast data has about 23 million nodes and 23 million edges, and building it takes several gigabytes of memory. Walking a copy of the graph would need memory for a second full copy, so by default assemble_contigs walks the original graph directly. This uses less memory, but it empties self.graph as the edges are used. This was pointed out during peer review: graph = self.graph doesn’t make a copy, so the original graph is changed. To make the trade-off explicit, assemble_contigs takes a copy_graph option. With copy_graph=True, it walks a copy and leaves self.graph unchanged. The final run peaked at 4.74 GB of memory on a 16 GB laptop.

**Different results from the same data.** Our runs gave very different results:
- 30 contigs with a longest contig of 22.6 million bp
- 7,632 contigs with a longest contig of 233,925 bp (Note: this is the version before our finalized version where we ended up removing reverse-complement duplicates.)
At first we didn't know why, or which was right. We found that our graphs were identical, because the total lengths differ by exactly 30 letters per extra contig, and that the difference came down to one word in the walk: `if` vs `while`.
With `while`, after reaching a dead end, the walk goes back to an earlier fork, takes the other branch, and adds it to the same path, even though the two branches are not connected.
With `if`, each walk is one straight path, and the other branch becomes its own contig. The 22.6 million bp contig is longer than the whole genome, which showed us that a Eulerian walk that uses every edge is not necessarily a realistic genome.

**The genome was assembled twice.** Our total assembly length was 192.1% of the genome size. We checked our data and found the reads come from both DNA strands about equally, so most regions were assembled twice: once as a contig and once as its reverse complement.
We added defremove_reverse_complement_duplicates(self, contigs, cutoff=0.9):going through the contigs longest first and dropping any contig where at least 90% of its k-mers are reverse complements of k-mers already kept. That brought the total down to 106.9% of the genome (2,967 contigs, N50 36,194 bp). The 90% cutoff is a judgment call. About 6.9% extra remains, which may come from shorter twins that differ by more than the cutoff allows, but we did not confirm this.

### Things I still don't understand
How real assemblers decide where to stop a contig at a fork, and how they handle reverse-strand duplicates while building the graph, instead of afterward like we did.

### Workflow
With busy schedules over the week, we struggled to find times when we were all available to meet. As a result, each of us individually thought of how to approach the problem and came together to compare our solutions. Our final submitted notebook is a combination of all the techniques we thought worked best from each individual approach. We discussed as a group where approaches differed and which sections to take from each, making sure that we were all on board with every decision. Although we weren't able to meet regularly throughout the week, we were able to make our collaboration style work for this project.

## Appendix (AI usage, citations, etc.)
We used Claude to help diagnose the Python recursion limit we ran into when running our notebooks against the real dataset for the first time. We prompted, "How can we increase the Python recursion limit? Will there be any unexpected consequences of this?" We used feedback from this prompt to increase the recursion limit to a suitable amount.
