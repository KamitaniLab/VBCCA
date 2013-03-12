% Variational Bayesian (VB) estimation for Bayesian Canonical Correlation Analysis
% This function is the main loop of VB estimation.
% 
% Generative model
%  x1  = W1*Z
%  x2  = W2*Z
%  x_i : observation (D_i x N), i = 1,2
%  W_i : weight (D_i x M), i = 1,2 
%  z   : latent variable (M x N)
%  b_i : bias (D_i x 1), i = 1,2
%  N   : number of samples
%  D_i : dimension of observation i
%  M   : dimension of latent variable
% Sparseness priors are assumed for all elements of W_i
% 
% [tr_struct,tr] = BCCAtrainMain(x1,x2,parm)
% 
% INPUT
% x1    : data 1 (D_1 x N)
% x2    : data 2 (D_2 x N)
%
% param : initial values of parameters (structure)
%  .M  = parm.M : initial dimensions of latent variable
%  .beta_inv{i} : observation noise i (1 x 1)
%  .A_inv{i}    : element-wise initial hyper-prior variance of W1 (D_i x M)
%  .A0_inv{i}   : element-wise hyper-prior A_inv{i} (D_i x M)
%  .gamma0{i}   : hyper-prior of confidence parameter gamma (D_i x M)
%  .thres_a_inv : threshold for a_inv (e.g., 1e-5) if 0, no thresholding
%  .Nitr        : iteration times
%  .NitrDisp    : how often to display iteration times
% 
% OUTPUT
% tr_struct     : structure of trained parameters
%  .W           : weight matrix (1x2 cell)
%  .W_inv       : projection matrix from x_i to z (1x2 cell)
%  .SigmaW_inv  : variances of W (1x2 cell)
%  .beta_inv    : inverse variance (1x2 cell)
%  .SigmaZ      : inverse covariance of z
%  .A_inv       : inverse variances of W (1x2 cell)
% 
% tr            : object including all data and estimated variables
%
function [tr_struct,tr] = BCCAtrainMain(x1,x2,parm)

if ~isfield(parm,'NitrDisp')
    parm.NitrDisp = 0;
end

% constructor (input observation set)
tr = BCCAtrain(x1,x2);
clear x1 x2

% initialize variables
initVar_common(tr,parm);
initVar_alpha_elementWise(tr,1,parm);
initVar_alpha_elementWise(tr,2,parm);

% VB iteration
for itr = 1:parm.Nitr
    if mod(itr,parm.NitrDisp)==0; disp(['Iteration: ' num2str(itr)]); end
    
    W_elementWise(tr,1);
    W_elementWise(tr,2);
    Zstep(tr);
    alpha_elementWise(tr,1);
    alpha_elementWise(tr,2);
    beta_step(tr,1);
    beta_step(tr,2);
end

tr_struct = save_parm(tr);