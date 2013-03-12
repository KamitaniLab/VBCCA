% Superclass of BCCA program
% 
% Methods
%  - obj = vbBCCA(Data)
%    Constructor. Allocate variables for input data.
%  This method is called by subclass functions BCCAtrain and BCCApred.
% 
classdef vbBCCA < handle
    properties
        x
        datDir
        datName
        mu
        Dx
        N
        M
        beta_inv
        W
        W_inv
        SigmaW_inv
        Z
        SigmaZ
        SigmaZ_inv
    end
    
    methods
        % constructor
        function obj = vbBCCA(Data)
    
            if length(Data)==2
                % load data as 2 input arguments                                
                obj.x{1} = Data{1};
                obj.x{2} = Data{2};
                
            else
                % load data from .mat directly
                ix_path = 0;
                for vn = 1:length(Data)
                    if strcmp(Data{vn},'dataPath')
                        ix_path    = vn+1;
                        obj.datDir = Data{ix_path};
                    end
                end
                
                if ix_path
                    for vn = 1:length(Data)
                        if strcmp(Data{vn},'x1')
                            load(Data{ix_path},Data{vn+1});
                            eval(['obj.x{1}=' Data{vn+1}]);
                            eval(['clear ' Data{vn+1}]);
                            obj.datName{1} = Data{vn+1};
                        end
                        if strcmp(Data{vn},'x2')
                            load(Data{ix_path},Data{vn+1});
                            eval(['obj.x{2}=' Data{vn+1}]);
                            eval(['clear ' Data{vn+1}]);
                            obj.datName{2} = Data{vn+1};
                        end
                    end
                    
                else
                    % load data as func('x1',x1,'x2',x2)
                    for vn = 1:length(Data)
                        if strcmp(Data{vn},'x1')
                            obj.x{1} = Data{vn+1};
                        end
                        if strcmp(Data{vn},'x2')
                            obj.x{2} = Data{vn+1};
                        end
                    end
                end
            end
            
            if ~isprop(obj,'x')
                error('Input must be (x1,x2) or ("dataPath","/home/...","x1",a,"x2",b)')
            end
            
            if ~isempty(obj.x{1})&&~isempty(obj.x{2})
                if size(obj.x{1},2)~=size(obj.x{2},2)
                    error('Sample size of x1 and x2 must be same');
                end
            end
            
            Ntmp = zeros(2,1);
            for nx = 1:length(obj.x)
                [obj.Dx{nx},Ntmp(nx)] = size(obj.x{nx});
            end
            obj.N = max(Ntmp);
            
        end
    end
end