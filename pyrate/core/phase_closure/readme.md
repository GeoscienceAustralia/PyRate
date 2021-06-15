## Phase Closure

To use the PyRate _phase closure_ correction, simply add the string _phase_closure_
to the list of `correct`ions as can be seen below:

```
[correct]
steps =
    orbfit
    refphase
    demerror
    phase_closure
    mst
    apscorrect
    maxvar
```

The following parameters are used from the PyRate config file:
```
# closure_thr: closure threshold for each pixel in multiples of pi, e.g. 0.5 = pi/2, 1 = pi.
# avg_ifg_err_thr: ifgs with more than this fraction of pixels above the closure threshold on average will be dropped entirely.
# min_loops_per_ifg: ifgs are dropped entirely if they do not participate in at least this many closure loops.
# max_loop_length: closure loops with up to this many edges will be used.
# max_loop_redundancy: A closure loop will be discarded if all constituent ifgs in that loop have
#                    already contributed to a number of loops equal to this parameter.
closure_thr:         0.5
avg_ifg_err_thr:     0.05
min_loops_per_ifg:   2
max_loop_length:     4
max_loop_redundancy: 2
```

The PyRate _phase closure_ correction has the following main functionalities:

1. Find the closed loops having edges less than or equal to _max_loop_length_.
   This is done in the function _mst_closure.sort_loops_based_on_weights_and_date_.
   We perform several steps in this stage:
    
    - Find closed loops with edges numbering between 3 and _max_loop_length_ (mst_closure.__find_closed_loops).
    - Sort ifgs within each loop based on first date (earlier date first). In case of a tie, we sort 
      based on the second date. We then compute weight of each ifg as the temporal baseline in days.
      Then we sum the interferogram weights in a loop to find the weight of each closure loop
      (mst_closure.__add_signs_and_weights_to_loops). 
    - Sort the loops based on weights, and in case of ties, further sort by primary date, and then by secondary
      date (mst_closure.sort_loops_based_on_weights_and_date).

2. Discard closure loops when all of the constituent ifgs have already contributed to _max_loop_redundancy_ loops.
   This is done in the function _closure_check.discard_loops_containing_max_ifg_count_.

3. Drop interferograms from further PyRate processing when they are found to not form part of any closed loop.
   This is done in the function _closure_check.__drop_ifgs_if_not_part_of_any_loop_.

4. Compute phase closure sums for each pixel (in radians) in all the chosen loops and flag when
   the resultant absolute sum exceeds _closure_thr_ * pi.
   Optionally the median closure sum across all pixels can be subtracted from the closure (parameter _subtract_median_; default to off)
   This is done in the function _sum_closure.__compute_ifgs_breach_count_.

5. Next, interferograms are dropped (removed from the processing list) if the fraction of constituent pixels
   breaching _closure_thr_ averaged over the loops the ifg participates in exceeds the _avg_ifg_err_thr_ 
   threshold, or the ifg does not contribute to a number of loops equal to or exceeding _min_loops_per_ifg_.
   This is done in the function _closure_check.__drop_ifgs_exceeding_threshold_.
   
6. Steps 1-5 are repeated iteratively until a stable list of interferograms are returned.
   The iteration is orchestrated by the function _closure_check.iterative_closure_check_.
   
7. Once a stable list of interferograms is found, in _correct.py_, write a new ifglist in the working directory, 
   update params, and use the updated ifglist for further PyRate processing.
   
8. Finally, find pixels in the interferogram phase data that breach _closure_thr_,
   and mask (assign NaNs) to those pixels in those interferograms.
   This is done in the function _closure_check.mask_pixels_with_unwrapping_errors_.

This completes the PyRate _phase closure_ correction.
