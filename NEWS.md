Release 1.1-r68 (4 September 2025)
----------------------------------

Notable changes:

 * Improvement: a theoretically better algorithm. On CHM13, the new algorithm
   finds 1448bp more LCRs than the old algorithm and misses 27bp out of 280Mb
   identified LCRs - the result is nearly the same. The performance is also
   similar.

 * Improvement: conform to the C99 standard.

 * Change: removed option -d for initial X-drop.

With X-drop disabled, the new algorithm is almost exact except a very rare
corner case where a full window is not LCR by definition but could be LCR if we
extend the window by a few basepairs. I don't know how to address this issue
without greatly reducing the performance. In practice, the issue only affects
<200bp (out of 280Mb identified LCRs) and will not impact downstream analysis
given the ambiguity in defining LCRs and the X-drop heuristic.

(1.1: 4 September 2025, r68)
