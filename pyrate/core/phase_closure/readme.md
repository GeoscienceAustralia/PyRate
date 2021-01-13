## Phase Closure

To use _phase closure_ correction, simply add the string _phase_closure_ in the list of `correct`ions as can be seen 
below:

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


The following params are used from the pyrate config file:
```
# phase closure params
# large_deviation_threshold_for_pixel # pi/2, pi=3.1416
# threshold_to_remove_ifg: ifgs with more than this fraction of pixels with error will be dropped
# loop_count_for_threshold_to_remove_ifg: pixel with phase unwrap error in at least this many loops
# phase_unwrap_error_threshold: pixel with phase unwrap error in more than this many ifgs will be flagged
# max_loop_length: loops upto this many edges are considered for closure checks
# subtract_median_in_closure_check: whether to subtract median during closure checks
# max_loop_count_for_each_ifgs: loops with more than these many ifgs are discarded.
# max_loop_count_for_each_ifgs: Ifg count must be met for all ifgs in the loop for loop to be discarded
large_dev_thr: 1.5708
threshold_to_remove_ifg: 0.07
loop_count_for_threshold_to_remove_ifg: 2
phase_unwrap_error_threshold: 5
max_loop_length: 3
subtract_median_in_closure_check: 1
max_loop_count_for_each_ifgs: 2
```

**Currently only max_loop_length: 3 produce repeatable results with Mexico CropA dataset.**

_Phase closure_ correction has the following main functionalities:

1. Compute the closed loops using `networkx`. Loops are assigned signs for each interferogram, and assigned a weight 
   based on total weight of each loop, which is the sum of difference between the ifg second and first date. This 
   is done in python file _mst_closure.py_. We perform several steps in this stage:
    
    1. Discard loops that are more than _max_loop_length_.
    2. Sort each loop based on first date of each interferogram (lower weight first). In case of a tie, we sort 
       based on the second date of the interferograms.
    3. Compute weight of each interferograms (=second date -first date). 
    3. Then we sum the weights of interferograms in a loop to find the weight of each closed loop.
    4. Sort the loops based on weights. In case of ties, we further sort by primary dates, and then by secondary 
       dates.
    5. Discard loops containing _max_loop_count_for_each_ifgs_. All ifgs in the loop must have contributed to at 
       least _max_loop_count_for_each_ifgs_. 
    6. Drop ifgs not part of any loop after this stage.

2. compute _sum_closure_ of each loop from stage 1 for each pixel. In addition, we perform the following steps in 
   _sum_closure.py_: 
    1. Find pixels breaching _large_dev_thr_. Note _large_dev_thr_ is specified in radians. Therefore at this stage 
       we need to convert phase data (in millimeters) into radians (check functions _shared.convert_to_radian_ and 
       it's use in the _Ifg_ class). 
    2. In order to find _Persistent Scatter_, compute the _check_ps_ for each pixel for each ifg.
    3. See use of _subtract_median_in_closure_check_ in function _sum_closure_.__compute_check_ps_.

3. _closure_check.py_ is used for orchestration of the functionalities above. After stage 2, we drop 
   ifgs exceeding _threshold_to_remove_ifg_ and _loop_count_for_threshold_to_remove_ifg_. See docstring in function 
   _closure_check.drop_ifgs_exceeding_threshold_.
   
4. Steps 1-3 are repeated until a stable list of interferograms are returned (see  
   _closure_check.filter_to_closure_checked_ifgs_).
   
5. Once a stable list of interferograms is found, in _correct.py_, write a new ifglist in the working directory, 
   update params, and use the updated ifglist for further pyrate processing.
   
6. Also in _correct.py_, we detect persistent scatterers (pixels) breaching _phase_unwrap_error_threshold_.

7. Finally, we write the ifg phase data after assigning nulls to pixels breaching _phase_unwrap_error_threshold_. 
   _Phase closure_ correction is done at this stage.