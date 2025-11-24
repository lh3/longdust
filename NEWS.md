Release 1.4-r97 (24 November 2025)
----------------------------------

Notable changes:

 * Improvement: almost doubled performance with tighter bound and 16-bit
   counts. The results are slightly different (2.5kb out of 277Mb).

(1.4: 24 November 2025, r97)



Release 1.3-r87 (16 November 2025)
----------------------------------

Notable changes:

 * New feature: correct for GC content. Biased GC only increases the scaling
   function f().

(1.3: 16 November 2025, r87)



Release 1.2-r75 (8 September 2025)
----------------------------------

Notable changes:

 * Bugfix: wrong formula. Correcting the formula slows down longdust by 60-80%,
   but errors like this need to be fixed. The output becomes ~1% smaller.

(1.2: 8 September 2025, r75)



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
