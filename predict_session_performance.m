function [r, p, pos_overlap, neg_overlap, predicted_behav] = predict_session_performance(fc_mats, behav, p_thresh, nrand)
% This function trains and tests a CPM with pre-calculated functional 
% connecitivty matrices using leave-one-out cross-validation. It calculates
% non-parametric significance when nrand>0.
% 
% INPUT
% fc_mats: Pre-calculated MxMxP matrix containing session-specific
% functional connectivity matrices, where M = number of network nodes and 
% P = number of sessions.
%
% behav: Px1 vector of scores for the behavior of interest
%
% p_thresh: p-value threshold for feature selection. In Rosenberg et al.,
% (2020), p_thresh = .01
%
% nrand: Number of randomizations for generating null distribution and
% determining non-parametric significance. Permutation testing is not run
% if rand=0. In Rosenberg et al. (2020), nrand = 1000.
%
% OUPTUT
% [r, p]: Output of Matlab corr function reflecting relationship between
% observed and predicted behavioral scores. If rand>0, p is a
% non-parametric p-value determined with permutation testing.
%
% pos_overlap: network mask of edges that predicted higher behavioral
% scores in every round of cross-validation
%
% neg_overlap: network mask of edges that predicted lower behavioral
% scores in every round of cross-validation
%
% predicted_behav: model predictions

% Set variables
[r, parametric_p, pos_overlap, neg_overlap, predicted_behav] = leave1out_cpm(fc_mats, behav, p_thresh);

if nrand>0
    rng(0,'twister');
    for i = 1:nrand
        [r_rand(i,1), ~, ~, ~, ~] = leave1out_cpm(fc_mats, behav(randperm(length(behav))), p_thresh);
    end
    p = (1+size(find(r_rand>=r),1))/(nrand+1); % get non-parametric p-value
else
    p = parametric_p;
end
