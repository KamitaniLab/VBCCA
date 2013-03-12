% Calculate Predictive Distributions for p(x2|x1) and p(x1|x2)
% 
% [pr_struct,pr] =BCCApredBoth(x1,x2,tr_struct,pr_parm)
% 
% INPUT
%  x1          : data 1
%  x2          : data 2
%  tr_struct   : structure of trained parameters by BCCAtrain
%
%  pr_parm     : input parameters to construct predictive distribution (structure)
%                .gvn_bias_flag{1 or 2} = 0: don't subtract bias of data
%                                         1: calculate bias of the test data and subtract it from the test data (1x2 cell) (default 1)
%                                         2: subtract bias of the training data from the test data
%                .SigmaZ_mode: 0: don't use the SigmaZ calculated in BCCAtrain
%                              1: use SigmaZ (default)
%  
% OUTPUT
%  pr_struct : structure of parameters of predictive distribution
%              .x_pr     : predicted value (mean of predictive distribution)
%              .x_pr_cov : covariance of predictive distribution
%              .prMat    : prediction matrix. x_j is predicted by equation x_j = prMat * x_i (where ix_gvn = i)
%              .bias     : bias of input data
% 
%  pr: object for calculating predictive distribution
%
function [pr_struct,pr] =BCCApredBoth(x1,x2,tr_struct,pr_parm)

pr = BCCApred('x1',x1,'x2',x2,'tr',tr_struct);    
clear x1 x2

% predict x1 from x2
pr_parm.ix_gvn = 2;
pred_simple(pr,pr_parm);

% predict x2 from x1
pr_parm.ix_gvn = 1;
pred_simple(pr,pr_parm);

% transform object pr to structure
pr_struct = save_parm(pr);


