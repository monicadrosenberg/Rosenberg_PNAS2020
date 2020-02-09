function [r, p, predicted_behav] = apply_sacpm(fc_mats, behav, path)
% This function applies the sustained attention connectome-based predictive
% model (saCPM) to data from novel individuals. The saCPM was defined in:
% Rosenberg MD, Finn ES, Scheinost D, Papademetris X, Shen X, 
% Constable RT & Chun MM. (2016). A neuromarker of sustained
% attention from whole-brain functional connectivity. Nature Neuroscience 
% 19(1), 165-171. Model coefficients were defined using robustfit in Matlab.
% 
% INPUT
% fc_mats: Pre-calculated MxMxP matrix containing individual-subject 
% functional connectivity matrices, where M = number of network nodes and 
% P = number of subjects.
%
% behav: Px1 vector of scores for the behavior of interest in the test set
% 
% path: path to saCPM.mat

% OUPTUT
% [r, p]: Output of Matlab corr function reflecting relationship between
% observed and predicted behavioral scores
%
% predicted_behav: saCPM predictions

% Load saCPM
load([path 'saCPM.mat'])

% Apply model to new data
for i = 1:size(fc_mats,3)
    net_strength.high(i,1) = sum(sum(high_attention_mask.*fc_mats(:,:,i)));
    net_strength.low(i,1)  = sum(sum(low_attention_mask.*fc_mats(:,:,i)));
    predicted_behav(i,1)   = robGLM_fit_update(1) + robGLM_fit_update(2)*(net_strength.high(i) - net_strength.low(i));
end

% Correlate predicted and observed behavior
[r, p] = corr(behav, predicted_behav, 'type','spearman');
