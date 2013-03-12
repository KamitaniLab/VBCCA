% This subclass offers METHODS to calculate predictive distributions
% 
% Method
% - pr = BCCApred(varargin)
%    Constructor. Allocate variables for data and parameters and make an object 'pr'
% - pred_simple(pr,parm)
%    Calculate a predictive distribution
% - pr = save_parm(pr)
%    Transform parameters of predictive distribution from the object 'pr' to struct 'pr' for saving
%
classdef BCCApred < vbBCCA
    properties
        x_pr
        x_pr_cov
        prZ
        prMat
        tr_mu;
        x_exist
    end
    
    methods
        % constructor
        function pr = BCCApred(varargin)
            % call superclass constructor
            pr = pr@vbBCCA(varargin);
            
            % initialize cells
            pr.mu       = cell(1,2);
            pr.prMat    = cell(1,2);
            pr.x_pr     = cell(1,2);
            pr.x_pr_cov = cell(1,2);
            pr.tr_mu    = cell(1,2);
            
            % find object of estimated parameters
            % load data from .mat directly
            for vn = 1:length(varargin)
                if strcmp(varargin{vn},'tr')
                    ix_tr = vn+1;
                end
            end

            pr.beta_inv{1}   = varargin{ix_tr}.beta_inv{1};
            pr.beta_inv{2}   = varargin{ix_tr}.beta_inv{2};
            pr.W{1}          = varargin{ix_tr}.W{1};
            pr.W{2}          = varargin{ix_tr}.W{2};
            pr.tr_mu{1}      = varargin{ix_tr}.mu{1};
            pr.tr_mu{2}      = varargin{ix_tr}.mu{2};
            pr.SigmaW_inv{1} = varargin{ix_tr}.SigmaW_inv{1};
            pr.SigmaW_inv{2} = varargin{ix_tr}.SigmaW_inv{2};
            pr.M             = varargin{ix_tr}.M;
            pr.Z             = [];
            pr.SigmaZ        = varargin{ix_tr}.SigmaZ;
            
            pr.x_exist       = find(~[isempty(pr.x{1}) isempty(pr.x{2})]);
        end
        
        % calculate predictive distribution
        function pred_simple(pr,parm)
            
            ix_gvn = parm.ix_gvn;
            if (ix_gvn~=1)&&(ix_gvn~=2)
                disp('ix_pr must be 1 or 2')
                return;
            end
            ix_pr = setxor([1 2],ix_gvn);
            
            % whether to calculate the bias of given data or not
            if isfield(parm,'gvn_bias_flag')
                if iscell(parm.gvn_bias_flag)
                    for n = pr.x_exist
                        if parm.gvn_bias_flag{n}==1
                            pr.mu{n} = mean(pr.x{n},2);
                            pr.x{n}  = pr.x{n} - repmat(pr.mu{n},[1 pr.N]);
                        elseif parm.gvn_bias_flag{n}==2
                            pr.x{n}  = pr.x{n} - repmat(pr.tr_mu{n},[1 pr.N]);
                        else
                            pr.mu{n} = [];
                        end
                    end
                else
                    if parm.gvn_bias_flag==1
                        pr.mu{ix_gvn} = mean(pr.x{ix_gvn},2);
                        pr.x{ix_gvn}  = pr.x{ix_gvn} - repmat(pr.mu{ix_gvn},[1 pr.N]);
                    elseif parm.gvn_bias_flag{n}==2
                        pr.x{ix_gvn}  = pr.x{ix_gvn} - repmat(pr.tr_mu{ix_gvn},[1 pr.N]);
                    else
                        pr.mu{ix_gvn} = [];
                    end
                end
            else % default setting
                pr.mu{ix_gvn} = mean(pr.x{ix_gvn},2);
                pr.x{ix_gvn}  = pr.x{ix_gvn} - repmat(pr.mu{ix_gvn},[1 pr.N]);
            end
            
            
            if isfield(parm,'SigmaZ_mode')&&(parm.SigmaZ_mode==0)
                SigmaZ = zeros(pr.M);
            else
                SigmaZ = pr.SigmaZ;
            end            
            
            %%% calculate posterior of z
            SigmaW2_inv_diag = diag(sum(pr.SigmaW_inv{ix_gvn},1));
            SigmaZnew        = SigmaZ + (1/pr.beta_inv{ix_gvn}).*(pr.W{ix_gvn}'*pr.W{ix_gvn} + SigmaW2_inv_diag);
            
            % covariance
            SigmaZnew_inv = inv(SigmaZnew);
            
            % prediction matrix
            pr.prZ{ix_pr}   = SigmaZnew_inv*((1/pr.beta_inv{ix_gvn}).*pr.W{ix_gvn}');
            pr.prMat{ix_pr} = pr.W{ix_pr}*pr.prZ{ix_pr};
            
            
            %%% predictive distribution
            % mean
            pr.x_pr{ix_pr} = pr.prMat{ix_pr}*pr.x{ix_gvn};
            % covariance
            pr.x_pr_cov{ix_pr} = pr.W{ix_pr}*SigmaZnew_inv*pr.W{ix_pr}' + pr.beta_inv{ix_pr}.*eye(size(pr.W{ix_pr},1));
        
            % whether to add bias to the predicted mean or not
            if isfield(parm,'pr_bias_flag')
                if iscell(parm.pr_bias_flag)
                    if parm.pr_bias_flag{ix_pr}==1
                        pr.x_pr{ix_pr} = pr.x_pr{ix_pr} + repmat(pr.tr_mu{ix_pr},[1 pr.N]);
                    end
                else
                    if parm.pr_bias_flag~=0
                        pr.x_pr{ix_pr} = pr.x_pr{ix_pr} + repmat(pr.tr_mu{ix_pr},[1 pr.N]);
                    end
                end
            end
        
        end
        %End pred_simple(pr,parm)
        
        
        % save parameters in structure
        function pr = save_parm(pr)
            pr_new.x = pr.x;
            pr_new.mu = pr.mu;
            pr_new.prZ = pr.prZ;
            pr_new.prMat = pr.prMat;
            pr_new.x_pr = pr.x_pr;
            pr_new.x_pr_cov = pr.x_pr_cov;
            pr = pr_new;
        end
    end
end