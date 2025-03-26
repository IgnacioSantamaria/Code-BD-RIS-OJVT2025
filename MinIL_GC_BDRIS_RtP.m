function [BDRIS_group,IL] = MinIL_GC_BDRIS_RtP(H,G,F,Mg,varargin)

% BD-RIS optimization for IL minimization
% The IL cost function does not take into account precoders/decoders
% We consider relaxed BD-RIS with power constraint trace(Theta'*Theta)<=M,
% which will be later projected into the set of unitary+symmetric matrices
%
% The symmetry can also be directly incorporated into the constraints
%
% Modification 17/02/25: iterative version of the group connected architecture


K = size(G,2);      % number of users
M = size(G{1},2);   % number of BDRIS elements

%% Default values
opt_params = struct();
opt_params.powerconstraint = 1;    % 1 -> tr(BDRIS'*BDRIS) <= M; 2-> globally passive (it depends on the precoders)
opt_params.symm = 1;               % 1-> symmetric BDRIS, 0-> unconstrained
opt_params.ini = 1;                % initialization identity matrix
opt_params.thresholdIL = 1e-20;
opt_params.MaxiterGC = 25;
if nargin < 4
    error(message('TooFewInputs'));
elseif nargin == 5
    params = varargin{1};
    for arg = fieldnames(params)'
        parameter = arg{1};
        param_value = params.(parameter);
        switch parameter
            case 'powerconstraint'
                opt_params.powerconstraint  = param_value;
            case 'symm'
                opt_params.symm  = param_value;
            case 'ini'
                opt_params.ini  = param_value;
            case 'thresholdIL'
                opt_params.thresholdIL  = param_value;
            case 'MaxiterGC'
                opt_params.MaxiterGC  = param_value;
        end
    end
elseif nargin > 5
    error(message('TooManyInputs'));
end

Groups = fix(M/Mg);        % number of groups

if opt_params.ini == 0
    BDRIS_group = zeros(M,M);  % we start with zeros
else
    BDRIS_group = eye(M);  % identity matrix (this is a feasible solution)
end
HBDRISaux = cell(K,K);
IL = 0;
niter = 0;
for i = 1:K  % tx
    for j = 1:K % rx
        HBDRISaux{i,j} = H{i,j} + F{j}*BDRIS_group*G{i}';  % Eq. channel for ith tx to jth rx
        if ne(i,j)
            IL =  IL + norm(HBDRISaux{i,j},'fro').^2;
        end
    end
end

true = 1;
while true == 1

    niter = niter+1;
    for gg = 1:Groups
        % make zero the group to update (to compute the fixed part of the
        % channel)
        BDRIS_group((gg-1)*Mg+1: gg*Mg,(gg-1)*Mg+1: gg*Mg) = zeros(Mg,Mg);
        % extract gth block matrix
        Fg = cell(1,K);
        Gg = cell(1,K);
        for kk=1:K
            Fg{kk} = F{kk}(:,(gg-1)*Mg+1: gg*Mg);  % block of Mg columns of F
            Gg{kk} = G{kk}(:,(gg-1)*Mg+1: gg*Mg);  % block of Mg columns of G
        end
        % Recompute the channel with the previous BDRIS block fixed
        Hg = cell(K,K);
        for i = 1:K  % tx
            for j = 1:K % rx
                Hg{i,j} = H{i,j} + F{j}*BDRIS_group*G{i}';  % Eq. channel for ith tx to jth rx
            end
        end
        % Compute gth block
        [BDRISg,~,~,~,~,~,~,~,~] = MinIL_BDRIS_closedform(Hg,Gg,Fg,opt_params);
        BDRIS_group((gg-1)*Mg+1: gg*Mg,(gg-1)*Mg+1: gg*Mg) = BDRISg;
    end
    HBDRISaux = cell(K,K);
    ILg = 0;
    for i = 1:K  % tx
        for j = 1:K % rx
            HBDRISaux{i,j} = H{i,j} + F{j}*BDRIS_group*G{i}';  % Eq. channel for ith tx to jth rx
            if ne(i,j)
                ILg =  ILg + norm(HBDRISaux{i,j},'fro').^2;
            end
        end
    end
    IL = [IL ILg];
    %% Check convergence
    if (abs(IL(end)-IL(end-1))<opt_params.thresholdIL)||(niter >= opt_params.MaxiterGC)
        true = 0;
    end
end