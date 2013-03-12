% This subclass offers METHODS to estimate parameters by variational Bayes (VB).
% It contains parameter initialization and estimation steps for VB
% 
% Methods
% - tr = BCCAtrain(varargin)
%    Constructor. Allocate variables for data and parameters and make an object 'tr'
% - initVar_common(tr,parm)
%    Initialize parameters (bias, Z, beta etc.)
% - initVar_alpha_elementWise(tr,ix,parm)
%    Initialize alpha parameters (A_inv, A0_inv)
% - initVar_alpha_single(tr,ix,parm)
%    Initialize alpha parameters (A_inv, A0_inv)
% - W_elementWise(tr,ix)
%    Estimate weight parameters (W_i) (in which a single alpha is assumed to a weight matrix W)
% - W_single(tr,ix)
%    Estimate weight parameters (W_i) (in which alphas are assumed to all elements of weight matrix W)
% - Zstep(tr)
%    Estimate latent variables (Z)
% - alpha_elementWise(tr,ix)
%    Estimate alpha parameters (A_inv) (in which a single alpha is assumed to a weight matrix W)
% - alpha_single(tr,ix)
%    Estimate alpha parameters (A_inv) (in which alphas are assumed to all elements of weight matrix W)
% - beta_step(tr,ix)
%    Estimate variance parameter beta
% - find_irrelative(tr)
%    Find irrelevant parameters of A_inv
% - tr = save_parm(tr)
%    Transform estimated parameters in the object 'tr' to struct 'tr' to save
%
classdef BCCAtrain < vbBCCA
    properties (SetAccess = private)
        Znorm_flag
        thres_a_inv
        SZZ
        SZZrep
        WW
        A_inv
        A0_inv
        gamma
        gamma0
        gamma_xx
        gamma_beta
        Meye
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function tr = BCCAtrain(varargin)
            % call superclass constructor
            tr            = tr@vbBCCA(varargin);
            
            % initialize superclass properties
            tr.beta_inv   = cell(1,2);
            tr.W          = cell(1,2);
            tr.W_inv      = cell(1,2);
            tr.SigmaW_inv = cell(1,2);
            tr.M          = 0;
            tr.Z          = 0;
            tr.SigmaZ     = 0;
            tr.SigmaZ_inv = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initialize variable
        function initVar_common(tr,parm)
            if nargin==1
                %%% default setting %%%
                % bias of observation
                for n = 1:2
                    tr.mu{n} = mean(tr.x{n},2);
                    tr.x{n}  = tr.x{n} - repmat(tr.mu{n},[1 tr.N]);
                end
                
                % choose smaller dimension
                if tr.Dx{1}<=tr.Dx{2}
                    tr.M = tr.Dx{1};
                else
                    tr.M = tr.Dx{2};
                end
                tr.beta_inv{1} = 1;
                tr.beta_inv{2} = 1;
                tr.thres_a_inv = 1e-7;
                tr.Znorm_flag  = 0;
                
            elseif nargin==2
                % whether to calculate bias or not
                if isfield(parm,'bias_flag')
                    for n = 1:2
                        if parm.bias_flag{n}
                            tr.mu{n} = mean(tr.x{n},2);
                            tr.x{n}  = tr.x{n} - repmat(tr.mu{n},[1 tr.N]);
                        else
                            tr.mu{n} = [];
                        end
                    end
                else % default setting
                    for n = 1:2
                        tr.mu{n} = mean(tr.x{n},2);
                        tr.x{n}  = tr.x{n} - repmat(tr.mu{n},[1 tr.N]);
                    end
                end
                
                tr.beta_inv{1} = parm.beta_inv{1};
                tr.beta_inv{2} = parm.beta_inv{2};
                tr.M           = parm.M;
                tr.thres_a_inv = parm.thres_a_inv;
                if isfield(parm,'Znorm_flag')
                    tr.Znorm_flag = parm.Znorm_flag;
                else
                    tr.Znorm_flag = 0; % default
                end
                
            else
                error('number of input must be 1 or 2')
            end
            
            % initialize Z
            tr.Z          = randn(tr.M,tr.N);
            tr.SigmaZ_inv = speye(tr.M);
            tr.SZZ        = tr.Z*tr.Z' + tr.N.*tr.SigmaZ_inv;
            tr.SZZrep{1}  = repmat(diag(tr.SZZ)',[tr.Dx{1} 1]);
            tr.SZZrep{2}  = repmat(diag(tr.SZZ)',[tr.Dx{2} 1]);
            
            % initialize other variables
            tr.Meye = speye(tr.M);
            tr.WW   = cell(1,2);
            for ix = 1:2
                tr.gamma_xx{ix}   = sum(sum(tr.x{ix}.^2))/2;
                tr.gamma_beta{ix} = tr.Dx{ix}*tr.N/2;
            end
        end
        
        function initVar_alpha_elementWise(tr,ix,parm)
            if nargin==2
                tr.A_inv{ix}  = ones(tr.Dx{ix},tr.M);
                tr.A0_inv{ix} = zeros(tr.Dx{ix},tr.M);
                tr.gamma0{ix} = zeros(tr.Dx{ix},tr.M);
                
            elseif nargin==3
                tr.A_inv{ix}  = parm.A_inv{ix};
                tr.A0_inv{ix} = parm.A0_inv{ix};
                tr.gamma0{ix} = parm.gamma0{ix};
                
            else
                error('number of inputs must be 1 or 2')
            end
            
            %initialize other variables
            tr.gamma{ix} = 1/2 + tr.gamma0{ix};
        end
        
        function initVar_alpha_single(tr,ix,parm)
            if nargin==2
                tr.A_inv{ix}  = 1;
                tr.A0_inv{ix} = 0;
                tr.gamma0{ix} = 0;
                
            elseif nargin==3
                tr.A_inv{ix}  = parm.A_inv{ix};
                tr.A0_inv{ix} = parm.A0_inv{ix};
                tr.gamma0{ix} = parm.gamma0{ix};
                
            else
                error('number of inputs must be 1 or 2')
            end
            
            %initialize other variables
            tr.gamma{ix} = tr.D{ix}*tr.M/2 + tr.gamma0{ix};
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % W-step
        function W_elementWise(tr,ix)
            tr.SigmaW_inv{ix} = tr.A_inv{ix}./((1/tr.beta_inv{ix}).*tr.SZZrep{ix}.*tr.A_inv{ix} + 1);
            tr.W{ix}          = (1/tr.beta_inv{ix}).*tr.x{ix}*tr.Z'.*tr.SigmaW_inv{ix};
            tr.WW{ix}         = diag(sum(tr.SigmaW_inv{ix},1)) + tr.W{ix}'*tr.W{ix};
        end
        
        function W_single(tr,ix)
            tr.SigmaW{ix}     = tr.beta{ix}.*tr.SZZ + ojb.Meye./tr.A_inv{ix};
            tr.SigmaW_inv{ix} = inv(tr.SigmaW{ix});
            tr.W{ix}          = tr.beta{ix}.*tr.x{ix}*tr.Z{ix}'*ojb.SigmaW_inv{ix};
            tr.WW{ix}         = tr.D{ix}.*tr.SigmaW_inv{ix} + tr.W{ix}'*tr.W{ix};
        end
        
        % Z-step
        function Zstep(tr)
            tr.SigmaZ     = (1/tr.beta_inv{1}).*tr.WW{1} + (1/tr.beta_inv{2}).*tr.WW{2} + tr.Meye;
            tr.SigmaZ_inv = inv(tr.SigmaZ);
            tr.W_inv{1}   = (1/tr.beta_inv{1}).*tr.SigmaZ_inv*tr.W{1}';
            tr.W_inv{2}   = (1/tr.beta_inv{2}).*tr.SigmaZ_inv*tr.W{2}';
            tr.Z          = tr.W_inv{1}*tr.x{1} + tr.W_inv{2}*tr.x{2};
            
            % normalization
            if tr.Znorm_flag
                Zstd = std(tr.Z,0,2);
                tr.Z = tr.Z./repmat(Zstd,[1 tr.N]);
            end
            
            % second moment
            tr.SZZ       = tr.Z*tr.Z' + tr.N.*tr.SigmaZ_inv;
            tr.SZZrep{1} = repmat(diag(tr.SZZ)',[tr.Dx{1} 1]);
            tr.SZZrep{2} = repmat(diag(tr.SZZ)',[tr.Dx{2} 1]);
        end
        
        % alpha-step
        function alpha_elementWise(tr,ix)
            tr.A_inv{ix} = (tr.W{ix}.^2./2 + tr.SigmaW_inv{ix}./2 + tr.gamma0{ix}.*tr.A0_inv{ix})./tr.gamma{ix};
        end
        
        function alpha_single(tr,ix)
            tr.A_inv{ix} = (sum(sum(tr.W{ix}.^2))./2 + tr.D{ix}.*sum(diag(SigmaW2_inv))./2 + tr.gamma0{ix}.*tr.A0_inv{ix})./tr.gamma{ix};
        end
        
        % beta-step
        function beta_step(tr,ix)
            gamma_zx        = trace(tr.W{ix}*tr.Z*tr.x{ix}');
            gamma_zzww      = trace(tr.SZZ*tr.WW{ix})/2;
            beta_inv_gamma  =  tr.gamma_xx{ix} - gamma_zx + gamma_zzww;
            tr.beta_inv{ix} = beta_inv_gamma./tr.gamma_beta{ix};
            %gamma1_zzww = trace((N.*SigmaZ_inv + Z*Z')*(diag(sum(SigmaW1_inv,1)) + W1'*W1))/2;
            %beta_inv_gamma1 = trace(x1*x1')/2 - trace(W1*Z*x1')+...
            %trace((N.*SigmaZ_inv + Z*Z')*(diag(sum(SigmaW1_inv,1)) + W1'*W1))/2;
        end
        
        % find irrelevance parameters
        function find_irrelative(tr)
            a1_inv     = sum(tr.A_inv{1},1);
            a2_inv     = sum(tr.A_inv{2},1);
            a1_inv_max = max(a1_inv);
            a2_inv_max = max(a2_inv);
            ix_a1      = find(a1_inv>a1_inv_max*tr.thres_a_inv);
            ix_a2      = find(a2_inv>a2_inv_max*tr.thres_a_inv);
            ix_z       = intersect(ix_a1,ix_a2);
            
            M0    = tr.M;
            tr.M  = length(ix_z);
            Mdiff = tr.M-M0;
            if Mdiff
                %disp([num2str(itr) ', cond of SigmaZ: ' num2str(condSZ) ': number of latent: ' num2str(M0) '->' num2str(M)]);
                %delete component
                tr.Z          = tr.Z(ix_z,:);
                tr.SigmaZ     = tr.SigmaZ(ix_z,ix_z);
                tr.SigmaZ_inv = inv(tr.SigmaZ);
                tr.Meye       = speye(tr.M);
                tr.SZZ        = tr.Z*tr.Z' + tr.N.*tr.SigmaZ_inv;
                tr.W{1}       = tr.W{1}(:,ix_z);
                tr.W{2}       = tr.W{2}(:,ix_z);
                tr.A_inv{1}   = tr.A_inv{1}(:,ix_z);
                tr.A_inv{2}   = tr.A_inv{2}(:,ix_z);
                tr.A0_inv{1}  = tr.A0_inv{1}(:,ix_z);
                tr.A0_inv{2}  = tr.A0_inv{2}(:,ix_z);
                tr.gamma{1}   = tr.gamma{1}(:,ix_z);
                tr.gamma{2}   = tr.gamma{2}(:,ix_z);
                tr.gamma0{1}  = tr.gamma0{1}(:,ix_z);
                tr.gamma0{2}  = tr.gamma0{2}(:,ix_z);
            end
        end
        
        % save parameters in structure
        function tr = save_parm(tr)
            tr_new.M          = tr.M;
            tr_new.W          = tr.W;
            tr_new.W_inv      = tr.W_inv;
            tr_new.SigmaW_inv = tr.SigmaW_inv;
            tr_new.beta_inv   = tr.beta_inv;
            tr_new.SigmaZ     = tr.SigmaZ;
            tr_new.A_inv      = tr.A_inv;
            tr_new.mu         = tr.mu;
            tr                = tr_new;
        end
    end
end