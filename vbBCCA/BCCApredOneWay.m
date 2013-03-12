% Calculate Predictive Distribution p(x2|x1) or p(x1|x2)
% 
% [pr_struct,pr] = BCCApredOneWay(x,tr_struct,pr_parm)
% 
% INPUT
%  x         : data (x1 of p(x2|x1) or x2 of p(x1|x2))
%  tr_struct : structure of trained parameters by BCCAtrain
%
%  pr_parm   : input parameters to construct predictive distribution (structure)
%              .ix_gvn : must be 1 or 2. If 1, p(x2|x1) is calculated. If 2, p(x1|x2) is calculated.
%              .gvn_bias_flag{1 or 2} = 0: don't subtract bias of data
%                                       1: calculate bias of the test data and subtract it from the test data (1x2 cell) (default 1)
%                                       2: subtract bias of the training data from the test data
%              .pr_bias_flag{i} = 0: don't add bias to predicted data (default) 
%                                 1: add bias of training data to predicted data
%              .SigmaZ_mode: 0: don't use SigmaZ calcuated in BCCAtrain
%                            1: use SigmaZ (default)
%  
% OUTPUT
%  pr_struct : structure of parameters of predictive distribution
%              .x_pr     : predicted value (mean of predictive distribution)
%              .x_pr_cov : covariance of predictive distribution
%              .prMat    : prediction matrix. x_j is predicted by equation x_j = prMat * x_i (where ix_gvn = i)
%              .bias     : bias of input data
% 
% pr: object for calculating predictive distribution
%
function [pr_struct,pr] = BCCApredOneWay(x,tr_struct,pr_parm)

if pr_parm.ix_gvn==1
    pr = BCCApred('x1',x,'x2',[],'tr',tr_struct);    
elseif pr_parm.ix_gvn==2
    pr = BCCApred('x1',[],'x2',x,'tr',tr_struct);    
end
clear x1 x2
pred_simple(pr,pr_parm);
pr_struct = save_parm(pr);
