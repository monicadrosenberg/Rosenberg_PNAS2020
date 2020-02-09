function [h, p, ci, stats, pct] = compare_network_strength(cond1_mats, cond2_mats, network_mask, n)
% This function calculates strength in a pre-calculated functional
% connectivity network in 2 conditions (e.g., rest and anesthesia). It
% compares network strength across conditions with a t-test and calculates
% the percent of the time this difference exceeds differences observed in 
% same-size random networks. 
% 
% INPUT
% cond1_mats: Pre-calculated MxMxP matrix containing all individual-
% subject functional connectivity matrices in condition 1, where M = number 
% of network nodes and P = number of subjects.

% cond2_mats: Pre-calculated MxMxP matrix containing all individual-
% subject functional connectivity matrices in condition 2, where M = number 
% of network nodes and P = number of subjects. Code assumes that the order 
% of subjects in cond1_mats matches the order in cond2_mats.

% network_mask: A pre-calculated MxM matrix containing 1s where a functional
% connection (i.e., edge) is included in the network and 0s elsewhere.
% Note: Code can be modified to calculate strength in multiple network
% simultaneously and compare to differences in network strength in
% non-overlapping same-size random networks. 

% n: Number of same-size random networks used to compare network strength
% difference. In Rosenberg et al., PNAS, n = 10000. If n = 0, comparison to
% random networks is not performed. 

% OUPTUT
% [h, pi, ci, stats]: Output of Matlab paired t-test function

% pct: Percent of random networks with smaller differences in network strength
% than the observed network.

% Calculate network strength for each subject and each condition
for s = 1:size(cond1_mats,3)
    net_strength.cond1(s,1) = nansum(nansum(cond1_mats(:,:,s).*network_mask));
    net_strength.cond2(s,1) = nansum(nansum(cond2_mats(:,:,s).*network_mask));
end

% Compare network strength across conditions
[h,p,ci,stats] = ttest(net_strength.cond1,net_strength.cond2);
                
% Compare network strength difference to same-size random networks
if n > 0 
    rng(0,'twister');
    node       = size(cond1_mats,1);
    aa         = ones(node,node);
    aa_upp     = triu(aa,1);
    upp_id     = find(aa_upp);  
    tstat_rand = zeros(n, 1);
    
    for randn = 1:n
        rand_net = zeros(node,node);
        tmp      = randi([1 (node*(node-1))/2],1,sum(sum(network_mask)));
        rand_net(upp_id(tmp(1:sum(sum(network_mask))/2))) = 1;
        rand_net = rand_net + rand_net';
        [~,~,~,stats_rand] = ttest(squeeze(nansum(nansum(rand_net.*cond1_mats))),squeeze(nansum(nansum(rand_net.*cond2_mats))));
        tstat_rand(randn)  = stats_rand.tstat;
    end
    pct = size(find(abs(tstat_rand)<abs(stats.tstat)),1)/n*100;
else
    pct = NaN;
end
