## Introduction

Longdust identifies long highly repetitive STRs, VNTRs, satellite DNA and other
low-complexity regions (LCRs) in a genome. It is motivated by and follows a
similar rationale to [SDUST][sdust]. Unlike SDUST which is limited to short
windows, longdust can find centromeric satellite and VNTRs with long repeat
units.

Longdust also overlaps with tandem repeat finders (e.g. [TRF][trf],
[TANTAN][tantan] and [ULTRA][ultra]) in functionality. Nonetheless, it is not
tuned for tandem repeats with two or three copies, but instead finds regions
without clear tandem structure. Longdust complements TRF etc to some extent.

## Comparison to other tools

Longdust finds 280.4Mb of LCRs from the T2T-CHM13 analysis set. 228.4Mb of them
overlap with satellites (plus ~5Mb flanking) annotated by the T2T consortium,
33.1Mb of the remainder (52.0Mb) overlap with TRF (2 7 7 80 10 50 500 -l12), and
14.7Mb of the rest (18.9Mb) with SDUST (-t30). Only 4.1Mb is left. Most
longdust LCRs are found by other tools collectively.

In comparison, TRF finds 244.0Mb of TRs with four or more copies. 98.1% of them
are identified by longdust as well. On the contrary, of 30.6Mb of TRs with less
than four copies, only 15.7% overlap with longdust LCRs. Longdust is not tuned
for TRs with low copy numbers by default. With 349Mb, TANTAN (-w500 -s.85)
finds the most TRs. 70.0Mb of them do not overlap with the union of T2T
satellite, TRF, SDUST and longdust. Even if we reduce the score threshold from
50 to 30 with TRF, 63.3Mb is still left. TANTAN seems to be finding distinct TRs.

On performance, longdust and TANTAN both ran for ~35 minutes. SDUST was the fastest
at 4 minutes only. They all used less than 1.5GB memory. TRF took nearly 13
hours for T2T-CHM13. Setting "2 5 7 80 10 30 2000 -l20" took 19 hours and 12GB;
reducing `-l` to 12 did not yield any output in 40 hours. TRF was much faster on
GRCh38 due to the lack of long satellite arrays.

## Algorithm

### Notations

Let $`\Sigma`$ be the DNA alphabet. $`x\in\Sigma^*`$ is a DNA string.
$`t\in\Sigma^k`$ is a $k$-mer. $`\kappa(x)\subset\Sigma^k`$ is the set of
$k$-mers in $x$. $`c_x(t)`$ is the number of $k$-mer $t$ in $x$.
Let $`|x|`$ be the length of $x$ and $\ell(x)=|x|-k+1$ the number of $k$-mers
in $x$.

Suppose we are working with one long genome string. We use closed interval
$`[j,i]`$ to represent the substring starting at $j$ and ending at $i$,
including the end points. We may use "interval" and "subsequence"
interchangeably.

### Defining LCRs

Given a sequence $x$, let $`S(x)\ge0`$ be its complexity score. $x$ is a
*perfect interval* if no substring of $x$ is scored higher than $`S(x)`$; $x$
is a *good interval* if no prefix or suffix of $x$ is scored higher.
SDUST and longdust differ in the formulation of $S(x)$ and the capability to
find all perfect/good intervals.

### The SDUST algorithm

SDUST scores the complexity of $x$ with
```math
S_D(x)=\frac{\sum_{t\in\kappa(x)}c_x(t)(c_x(t)-1)/2}{\ell(x)}
```
It hardcodes $k=3$ and finds *all* perfect intervals of length $`\le`$64bp and
score $`\ge T`$ in a genome. It is an *exact* algorithm.

A sequence of length $w$ may have $O(w^2)$ perfect intervals in the worst case.
SDUST may traverse these intervals when processing a position. The overall time
complexity is thus $`O(w^3L)`$ where $w$ is the window length and $L$ is the
genome size. The cubic factor makes SDUST impractical for large $w$.

### The longdust algorithm

Longdust scores complexity with
```math
S_L(x)=\sum_{t\in\kappa(x)}\log\,c_x(t)!-f\left(\frac{\ell(x)}{4^k}\right)
```
where
```math
f(\lambda)=e^{-\lambda}\sum_{n=0}^\infty\log(n!)\cdot\frac{\lambda^n}{n!}
```
is calculated numerically. The first term in $`S_L(x)`$ comes from the log
probability of $x$ under a composite likelihood model. Please see the [math
notes](tex/notes.tex) for derivation. Given a threshold $`T\gt0`$, introduce
```math
S'_L(x,T)=S_L(x)-T\cdot\ell(x)
```
Longdust identifies a good interval $`[j,i]`$ via a backward and then a forward
scan through $`[i-w,i]`$ at each genomic position $i$. The time complexity is
$`O(wL)`$.  This algorithm only finds *a subset* of good intervals under
$`S'_L(x)`$ and is thus *heuristic*.

In the code, longdust impements a few strategies to speed up the search without
changing the output. It also uses BLAST-like X-drop to break at long non-LCR
intervals. Due to heuristics, longdust may generate slightly different output
on the reverse complement of the input sequence. For strand symmetry like
SDUST, longdust takes the union of intervals identified from both strands.

[sdust]: https://pubmed.ncbi.nlm.nih.gov/16796549
[trf]: https://github.com/Benson-Genomics-Lab/TRF
[ultra]: https://github.com/TravisWheelerLab/ULTRA
[tantan]: https://gitlab.com/mcfrith/tantan
