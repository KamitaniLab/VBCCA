% This script demonstrates parameter estimation and visual image reconstruction by Bayesian CCA.
% Model parameters are trained by random visual images and corresponding fMRI activity patterns.
% After that, visual images of letters and geometric shapes are reconstructed from fMRI data.

addpath('../vbBCCA/');
train_file = '../data/V1_raw_random.mat';
test_file = '../data/V1_mean_figure.mat';

load(train_file,'I','R');

%%% parameter settings for Bayesian CCA
D1                  = size(I,1);
D2                  = size(R,1);
M                   = 100; %dimension of latent variable z
tr_parm.M           = M;
tr_parm.beta_inv{1} = 1;
tr_parm.beta_inv{2} = 1;
tr_parm.A_inv{1}    = ones(D1,M);
tr_parm.A_inv{2}    = ones(D2,M);
tr_parm.A0_inv{1}   = zeros(D1,M);
tr_parm.A0_inv{2}   = zeros(D2,M);
tr_parm.gamma0{1}   = zeros(D1,M);
tr_parm.gamma0{2}   = zeros(D2,M);
tr_parm.thres_a_inv = 0;%1e-6;
tr_parm.Nitr        = 1000;
tr_parm.NitrDisp    = 100;
tr_struct           = BCCAtrainMain(I,R,tr_parm);

%%% plot estimated weight matrix
a1_inv              = sum(tr_struct.A_inv{1},1);
[tmp,ix]            = sort(a1_inv,'descend');

W = tr_struct.W{1};

NR = 10;
NC = 10;
figure;
for n = 1:NR*NC;
    subplot(NR,NC,n)
    imagesc(reshape(W(:,ix(n)),10,10));
    setfigure;
    caxis([-0.5 0.5]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% predict images from brain activity  %%%
load(test_file,'I2','R2');
I = I2; R = R2;
pr_parm.SigmaZ_mode = 0;
pr_struct = BCCApredBoth(I,R,tr_struct,pr_parm);
Ip = pr_struct.x_pr{1};

% normalized by sigmoid function
Ipn = 1./(1+ exp(-10.*Ip));

% calculate reconstruction error
Error = mean((Ipn-I).^2,1);

%%% plot reconstructed images %%%
NC1 = 5;
NR1 = 8 +1;

%%% figures
ix_tmp1(1,:) = [1 9 11 18 21 29 31 40];  %square
ix_tmp1(2,:) = [2 10 12 16 22 30 32 37]; %square open
ix_tmp1(3,:) = [3 6 13 20 23 28 33 39];  %cross
ix_tmp1(4,:) = [4 8 14 17 24 27 34 38];  %X
ix_tmp1(5,:) = [5 7 15 19 25 26 35 36];  %square large

% ascending order for reconstruction error for each figure
for n = 1:size(ix_tmp1,1)
    [e_tmp, ix_tmp0] = sort(Error(ix_tmp1(n,:)));
    ix_tmp1(n,:) = ix_tmp1(n,ix_tmp0);
end
clear ix_tmp0
ix1 = ix_tmp1(:);

figure;
% presented    
for n = 1:5;
    subplot(NR1,NC1,n);    
    imagesc(reshape(I(:,n),10,10));
    setfigure;
    caxis([0 1]);    
end

% predicted
for n = 1:length(ix1);
    subplot(NR1,NC1,n + NC1)
    imagesc(reshape(Ipn(:,ix1(n)),10,10));
    setfigure;    
    caxis([0 1]);    
end
S= get(0,'ScreenSize');
Hight = round(S(4)*0.8);
Width = Hight*2/3;
Left1 = round(S(3)/2-Width);
Bottom = round(S(4)*0.1);
set(gcf,'Position',[Left1 Bottom Width Hight]);


%%% characters

% ascending order for reconstruction error for each character
% thin characters
ix_tmp2(1,:) = [41 49 51 58 61 69 71 80]; %n
ix_tmp2(2,:) = [40 42 52 56 62 70 72 77]; %e
ix_tmp2(3,:) = [43 46 53 60 63 68 73 79]; %u
ix_tmp2(4,:) = [44 48 54 57 64 67 74 78]; %r
ix_tmp2(5,:) = [45 47 55 59 65 66 75 76]; %o

for n = 1:size(ix_tmp2,1)
    [e_tmp, ix_tmp0] = sort(Error(ix_tmp2(n,:)));
    ix_tmp2(n,:) = ix_tmp2(n,ix_tmp0);
end
ix2 = ix_tmp2(:);

figure;
% presented    
for n = 1:5;
    subplot(NR1,NC1,n);    
    imagesc(reshape(I(:,40+n),10,10));    
    setfigure;        
    caxis([0 1]);    
end

% predicted
for n = 1:length(ix2);
    subplot(NR1,NC1,n + NC1)
    imagesc(reshape(Ipn(:,ix2(n)),10,10));
    setfigure;
    caxis([0 1]);
end
Left2 = round(S(3)/2);
set(gcf,'Position',[Left2 Bottom Width Hight]);
