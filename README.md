# Rosenberg_PNAS2020
Code associated with Rosenberg et al., 2020, PNAS

The included MATLAB scripts were used in the study: Rosenberg et al., 2020. Functional connectivity predicts changes in attention observed across minutes, days, and months. PNAS. https://www.pnas.org/content/early/2020/02/03/1912226117

(1) apply_saCPM.m: Applies the sustained attention connectome-based predictive model (saCPM) to novel participants' functional connectivity matrices to predict sustained attention function. Loads saCPM.mat. Used in Rosenberg et al. (2020) Experiments 2 & 3.

(2) predict_taskblock_performance.m: Trains a CPM using data from all but one participant, and applies it to multiple functional connectivity matrices from the left-out individual to predict multiple measures of behavior. Used in Rosenberg et al. (2020) Experiment 1.

(3) predict_session_performance.m: Trains and tests a CPM with pre-calculated functional connectivity matrices using leave-one-out cross validation and calculates non-parametric significance. Calls leave1out_cpm.m. Used in Rosenberg et al. (2020) Experiment 3, "Within-subject model comparison" section.

(4) compare_network_strength.m: Measures strength in a pre-calculated functional connectivity network in 2 conditions (e.g., rest and anesthesia) and compares network strength across conditions. Used in Rosenberg et al. (2020) Experiments 4 & 5.
