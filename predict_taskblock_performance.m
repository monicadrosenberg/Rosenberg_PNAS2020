function [r, p, null_r, pos_overlap, neg_overlap, predicted_behav] = predict_taskblock_performance(train_mats, test_mats, behav_train, behav_test, p_thresh, nrand)
% This function trains and tests a connectome-based predictive model (CPM)
% with pre-calculated functional connecitivty matrices using leave-one-out 
% cross-validation. The CPM is trained using data from all but one 
% participant, and applied to multiple functional connectivity matrices 
% from the left-out individual to predict multiple measures of behavior. 
% The function calculates non-parametric significance when nrand>0.
% 
% See Finn et al. (2015) Nature Neuroscience & Shen et al. (2017) Nature 
% Protocols for more detail on CPM. A CPM code repository is updated 
% and maintained at https://github.com/YaleMRRC/CPM.
%
% INPUT
% train_mats: Pre-calculated MxMxP matrix containing functional connectivity
% matrices, where M = number of network nodes and P = number of subjects.  
%
% test_mats: Pre-calculated Px1 cell where P = number of subjects. Each cell 
% contains a MxMxQ functional matrix, where M = number of network nodes and 
% Q = number of observations (i.e., task blocks) per subject. 
%
% behav_train: Px1 vector of scores for the behavior of interest (i.e.,
% overall task performance).
%
% behav_test: Px1 cell. Each cell includes a Qx1 vector of one subject's 
% scores for the behavior of interest (i.e., each subject's task block 
% performance).
%
% p_thresh: p-value threshold for feature selection. In Rosenberg et al.,
% (2020), p_thresh = .01
%
% nrand: Number of randomizations for generating null distribution and
% determining non-parametric significance. Permutation testing is not run
% if rand=0. In Rosenberg et al. (2020), nrand = 1000.
%
% OUPTUT
% [r, p]: Output of Matlab corr function reflecting within-subject 
% relationships between observed and predicted behavioral scores. If rand>0,
% p is determined with permutation testing.
%
% pos_overlap: network mask of edges that predicted higher behavioral
% scores in every round of cross-validation
%
% neg_overlap: network mask of edges that predicted lower behavioral
% scores in every round of cross-validation
%
% predicted_behav: model predictions for each task block and each subject

% Set variables
rng(0,'twister');
nmat   = size(train_mats,3); % number of matrices
node   = size(train_mats,1); % number of nodes
aa     = ones(node,node);
aa_upp = triu(aa,1);
upp_id = find(aa_upp);    % indices of edges in the upper triangular of an node x node matrix
n_edge = length(upp_id);  % total number of edges

pos_mask_all = zeros(node, node, nmat);
neg_mask_all = zeros(node, node, nmat);

for i = 1:nmat
    i
    
    % exclude data from left-out session
    train_mats_tmp        = train_mats;
    train_mats_tmp(:,:,i) = [];
    train_behav           = behav_train;
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
    for j = 1:size(test_mats{i},3)
        tmp_mat = test_mats{i};
        test_pos_mean{i,1}(j,1)   = nansum(nansum(pos_mask.*test_mats{i}(:,:,j)));
        test_neg_mean{i,1}(j,1)   = nansum(nansum(neg_mask.*test_mats{i}(:,:,j)));
        predicted_behav{i,1}(j,1) = b_train(1) + (b_train(2)*(test_pos_mean{i}(j) - test_neg_mean{i}(j)));  
    end
    
    [r(i,1), parametric_p(i,1)] = corr(behav_test{i}, predicted_behav{i}, 'type','spearman');
    
    if nrand>0
        for k = 1:nrand
            behav_rand = predicted_behav{i};
            behav_rand = behav_rand(randperm(length(behav_rand)));
            [null_r(i,k), ~] = corr(behav_rand, predicted_behav{i},'type','spearman');
        end
        p(i,1) = (1+size(find(null_r(i,:)>=r(i)),2))/(nrand+1);
    else
        p = parametric_p;
    end
end

pos_overlap = zeros(node,node);
pos_overlap(sum(pos_mask_all,3)==nmat) = 1;
neg_overlap = zeros(node,node);
neg_overlap(sum(neg_mask_all,3)==nmat) = 1;
