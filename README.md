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

### The SDUST algorithm

Given a sequence $x$, SDUST defines its complexity score with
```math
S(x)=\frac{\sum_{t\kappa(x)}c_x(t)(c_x(t)-1)/2}{\ell(x)}
```
where $`\kappa(x)`$ is the set of $k$-mers in $x$, $`c_x(t)`$ is the number
of $k$-mer $t$ in $x$ and $`\ell(x)=|x|-k+1`$ is the number of $k$-mers. SDUST
by default uses $k=3$ and $|x|\le64$. $x$ is a *perfect interval* if the score
of any subsequence is no greater than $S(x)$.

[sdust]: https://pubmed.ncbi.nlm.nih.gov/16796549
[trf]: https://github.com/Benson-Genomics-Lab/TRF
[ultra]: https://github.com/TravisWheelerLab/ULTRA
