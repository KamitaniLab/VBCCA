% This script demonstrates parameter estimation and brain activity prediction by Bayesian CCA.
% Ten-fold cross validation is conducted using random images and corresponding fMRI activity patterns.
% Visual images are identified by comparing the predicted and measured fMRI data.

addpath('../vbBCCA/');
data_file = '../data/V1_raw_random.mat';

load(data_file,'I','R');

%%% parameter settings for training
D1                   = size(I,1);
D2                   = size(R,1);
M                    = 100; %dimension of latent variable z
tr_parm.M            = M;
tr_parm.beta_inv{1}  = 1;
tr_parm.beta_inv{2}  = 1;
tr_parm.A_inv{1}     = ones(D1,M);
tr_parm.A_inv{2}     = ones(D2,M);
tr_parm.A0_inv{1}    = zeros(D1,M);
tr_parm.A0_inv{2}    = zeros(D2,M);
tr_parm.gamma0{1}    = zeros(D1,M);
tr_parm.gamma0{2}    = zeros(D2,M);
tr_parm.thres_a_inv  = 0;
tr_parm.Nitr         = 1000;
tr_parm.NitrDisp     = 100;

%%% parameter settings for predictive distribution
pr_parm.ix_gvn       = 1; % to predict fMRI data from visual images
pr_parm.pr_bias_flag = 1; % add bias to predicted data, to compare predicted and measured data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10-fold cross validataion
 Ncv         = 1;
 Nblock      = 132;
 ixAll       = 1:1320;
 ixSet       = [1:Nblock:1320-Nblock+1;
                   Nblock:Nblock:1320];
 NpseudImage = 999;
 Ncandidate  = NpseudImage + 1;

C = [];
for nc = 1:Ncv;
    disp(['Fold: ' num2str(nc)])
    
    ix_test  = ixSet(1,nc):ixSet(2,nc);
    ix_train = setxor(ixAll,ix_test);
    
    % training data set
    Itrain = I(:,ix_train);    
    Rtrain = R(:,ix_train);        
    
    % test data set
    Itest = I(:,ix_test);
    Rtest = R(:,ix_test);
    
    % training
    tr_struct = BCCAtrainMain(Itrain,Rtrain,tr_parm);
    %Rtest = Rtest -repmat(tr_struct.mu{2},[1 Nblock]);
    
    % test (predict brain activities from presented visual images)
    pr_struct = BCCApredOneWay(Itest,tr_struct,pr_parm);
    
    Ctmp2 = zeros(Nblock,Ncandidate);
    for nb = 1:Nblock
        % generate pseudo visual images
        Ipseud = round(rand(D1,NpseudImage));
        
        % predict brain activities from pseudo visual images
        pr_struct_pseud = BCCApredOneWay(Ipseud,tr_struct,pr_parm);
                
        R_candidate = [pr_struct.x_pr{2}(:,nb) pr_struct_pseud.x_pr{2}];
        
        % calculate correlation between predicted and measured brain activity patterns
        for np = 1:Ncandidate
            Ctmp1 = corrcoef(Rtest(:,nb),R_candidate(:,np));
            Ctmp2(nb,np) = Ctmp1(1,2);
        end
    end
    C = cat(1,C,Ctmp2);
end

%calculate percent correct
Nstart = 10;
for Nend = Nstart:Ncandidate;
    [tmp,ix_max] = max(C(:,1:Nend),[],2);
    Ncorrect(Nend-Nstart+1) = length(find(ix_max==1));
end
Pcorrect = Ncorrect./size(C,1);

figure;
plot(Pcorrect)
xlim([1 Ncandidate-Nstart+1])
set(gca,'Xtick',[1 191 391 591 791 991]);
set(gca,'XTickLabel',{'10','200','400','600','800','1000'})
xlabel('Set Size')
ylabel('% correct')
set(gcf,'color','white');
set(gca,'Box','off');
set(gca,'TickDir','out');





