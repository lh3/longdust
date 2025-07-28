## Introduction

Longdust identifies long STRs, VNTRs, satellite DNA and other low-complexity
regions (LCRs) in a genome. In comparison to [SDUST][sdust] which use 3-mers in
64bp windows to find short LCRs, longdust by default uses 7-mers in 5kb windows
and can find centromeric repeats and VNTRs with long repeat units. It uses a
different algorithm as SDUST is impractical given long windows.

Longdust also overlaps with tandem repeat finders (e.g. [TRF][trf] and
[ULTRA][ultra]) in functionality. Nonetheless, it misses tandem repeats with
two or three copies but finds regions without clear tandem structure. It also
identifies high-order repeats (HORs) with kb-long repeat units.

In general, longdust does not replace SDUST or tandem repeat finders. It has a
different focus that other methods are not optmized for.

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

### The SDUST algorithm

Given a sequence $x$, SDUST scores its complexity with
```math
S_D(x)=\frac{\sum_{t\in\kappa(x)}c_x(t)(c_x(t)-1)/2}{\ell(x)}
```
It defines $x$ as a *perfect interval* if the score of any subsequence is no
greater than $`S_D(x)`$. SDUST hardcodes $k=3$ and reports all perfect
intervals of length $`\le`$64bp and score $`\ge T`$ in a genome.

In a sequence of length $w$, there are $O(w^2)$ perfect intervals in the worst
case. SDUST may traverse these intervals when processing a position. The overall
time complexity is thus $`O(w^3L)`$ where $w$ is the window length and $L$ is
the genome size. The cubic factor makes SDUST impractical for large $w$.

### The longdust algorithm

Longdust scores complexity with
```math
S_L(x)=\sum_{t\in\kappa(x)}\log\,c_x(t)!-f\left(\frac{\ell(x)}{4^k}\right)
```
where
```math
f(\lambda)=e^{-\lambda}\sum_{n=0}^\infty\log(n!)\cdot\frac{\lambda^n}{n!}
```
is calculated numerically. Please see the [math notes](tex/notes.tex) for
derivation. Given a threshold $`T\gt0`$, introduce
```math
S'_L(x,T)=S_L(x)-T\cdot\ell(x)
```
Longdust reports $`[j,i]`$ such that $`S'_L([j,i],T)>0`$ and no suffix or
prefix of $`[j,i]`$ is scored higher. This is different from SDUST which also
requires no internal substring can be scored higher.

Longdust finds $`[j,i]`$ via a backward and then a forward scan through
$`[i-w,i]`$ at each genomic position $i$. The time complexity is $`O(wL)`$.
It additionally implements a few strategies to speed up the search without
changing the output. It also uses BLAST-like X-drop to break at long non-LCR
intervals. The longdust algorithm would generate slightly different output on
the reverse complement of the input sequence. For strand symmetry like SDUST,
longdust takes the union of intervals identified from both strands.

[sdust]: https://pubmed.ncbi.nlm.nih.gov/16796549
[trf]: https://github.com/Benson-Genomics-Lab/TRF
[ultra]: https://github.com/TravisWheelerLab/ULTRA
