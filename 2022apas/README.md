# Workflow

## Find list of sessions for each stage
* `generate_sessions_each_mouse.py`: (~14s) Saves figuresdata/2022apas/sessions/pamoXXX.csv for each mouse.
  * Currently, we need to run separately for each cohort (2, 3 and 4)

## Calculate performance for each session/mouse for stage 3:
* `generate_learning_curve_stage3.py`: (~6s) Saves figuresdata/2022apas/learning_curve_stage3/fraction_correct_stage3.csv

## Figure 1
* `generate_example_performance.py`: save learning and psychometric for a single mouse.
* `figure_example_performance.py`

## Linear fit to learning
* studyutils.fit_learning_curves(dframe) # dframe comes from fraction_correct_stage3.csv

## Figure 2
* `generate_trials_per_stage.py`: total trials and trials per session for each stage.
* `generate_distributions_learning.py`: calculate distributions of P21 and T70
* `generate_learning_comparison.py`: calculate learning measurements and average learning curves.
* `figure_learning_comparison.py`
### To generate results for slow learners
* Re-run `generate_learning_comparison.py` with `learnerGroup = 'slow'`
### To generate results including newer mice
* Re-run `generate_distributions_learning.py` including cohorts 2,3,4.
* Re-run `generate_learning_comparison.py` using aoc=[2,3,4] in `studyutils.mice_each_condition()`

## Figure 3
* `generate_psychometric_each_mouse.py`: (~2sec)
run -t generate_psychometric_each_mouse.py early 2
run -t generate_psychometric_each_mouse.py early 3
run -t generate_psychometric_each_mouse.py late 2
run -t generate_psychometric_each_mouse.py late 3


## OBSOLETE: Calculate performance for each session/mouse for stage 4:
* `generate_learning_curve_stage4.py`: (~81s) Saves figuresdata/2022apas/learning_curve_stage4/fraction_correct_stage3.csv

# Extras
* `figure_learning_curve_each_mouse.py`: Plot learning curve for each mouse (excluding antibias sessions)
* ``

TO DO:
* DONE figure_example_performance: change 'fractionHitsEachValue' to 'fractionLeftEachValue'
* DONE Same for whatever script creates these data. example_performance_{}.npz'
* Rename generate_learning_comparison to old?  and other script?

