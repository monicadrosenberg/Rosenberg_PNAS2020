function [r, parametric_p, pos_overlap, neg_overlap, predicted_behav] = leave1out_cpm(fc_mats, behav, p_thresh)
% Trains and tests a connectome-based predictive model using leave-one-out
% cross-validation. See Shen et al. (2017) Nature Protocols for more detail.
% A CPM code repository is updated and maintained at
% https://github.com/YaleMRRC/CPM.
% 
% INPUT
% fc_mats: Pre-calculated MxMxP matrix containing session-specific
% functional connectivity matrices, where M = number of network nodes and 
% P = number of sessions.
%
% behav: Px1 vector of scores for the behavior of interest
%
% p_thresh: p-value threshold for feature selection. In Rosenberg et al.,
% 2020, PNAS, p_thresh = .01
%
% OUPTUT
% [r, parametric_p]: Output of Matlab corr function reflecting relationship 
% between observed and predicted behavioral scores
%
% pos_overlap: network mask of edges that predicted higher behavioral
% scores in every round of cross-validation
%
% neg_overlap: network mask of edges that predicted lower behavioral
% scores in every round of cross-validation
%
% predicted_behav: model predictions

% Set variables
nmat   = size(fc_mats,3); % number of matrices
node   = size(fc_mats,1); % number of nodes
aa     = ones(node,node);
aa_upp = triu(aa,1);
upp_id = find(aa_upp);    % indices of edges in the upper triangular of an node x node matrix
n_edge = length(upp_id);  % total number of edges

pos_mask_all = zeros(node, node, nmat);
neg_mask_all = zeros(node, node, nmat);

for i = 1:nmat
    i
    
    % exclude data from left-out session
    train_mats_tmp        = fc_mats;
    train_mats_tmp(:,:,i) = [];
    train_behav           = behav;
    train_behav(i)        = [];
    
    % create n_train_sub x n_edge matrix
    train_vect = reshape(train_mats_tmp, node*node, nmat-1)';
    upp_vect   = train_vect(:,upp_id);
    
    % relate behavior to edge strength in training data
    cp = zeros(n_edge, 1);
    cr = zeros(n_edge, 1);
    
    for ii = 1:n_edge
        [b,stats] = robustfit(upp_vect(:,ii), train_behav);
        cp(ii)    = stats.p(2);
        cr(ii)    = sign(stats.t(2))*sqrt((stats.t(2)^2/(nmat-1-2))/(1+(stats.t(2)^2/(nmat-1-2))));
    end
    
    % select edges based on threshold
    pos_edge = zeros(1, n_edge);
    neg_edge = zeros(1, n_edge);
    
    cp_pos           = find(cp<p_thresh & cr>0);
    pos_edge(cp_pos) = 1;
    cp_neg           = find(cp<p_thresh & cr<0);
    neg_edge(cp_neg) = 1;
    
    pos_mask = zeros(node, node);
    neg_mask = zeros(node, node);
    
    pos_mask(upp_id) = pos_edge;
    neg_mask(upp_id) = neg_edge;
    
    % save masks from every iteration of the leave-one-out loop
    pos_mask_all(:,:,i) = pos_mask;
    neg_mask_all(:,:,i) = neg_mask;
    
    for k = 1:nmat-1
        train_pos_mean(k,1) = nansum(nansum(pos_mask.*train_mats_tmp(:,:,k)));
        train_neg_mean(k,1) = nansum(nansum(neg_mask.*train_mats_tmp(:,:,k)));
    end
    
    % build model with training data
    b_train = robustfit(train_pos_mean-train_neg_mean,train_behav);
    
    % generate predictions for test data
    test_pos_mean(i,1) = nansum(nansum(pos_mask.*fc_mats(:,:,i)));
    test_neg_mean(i,1) = nansum(nansum(neg_mask.*fc_mats(:,:,i)));
    predicted_behav(i,1) = b_train(1) + (b_train(2)*(test_pos_mean(i) - test_neg_mean(i)));  
end

% Find egdes appearing in every round of leave-one-out cross-validation
pos_overlap = zeros(node,node);
neg_overlap = zeros(node,node);
pos_overlap(sum(pos_mask_all,3)==nmat) = 1;
neg_overlap(sum(neg_mask_all,3)==nmat) = 1;

% correlate predicted and observed behavior
[r, parametric_p] = corr(behav, predicted_behav,'type','spearman');