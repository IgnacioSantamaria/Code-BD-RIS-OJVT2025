function [BDRIS,IL, ILconv] = MinIL_GC_BDRIS(H,G,F,Mg,varargin)

% Description: Min-IL optimization of a Group-Connected BD-RIS architecture
% This is Algorithm 3 in [1]. For optimizing each group we can use either
% the MO or the (suboptimal) RtP algorithm.
%
% Input parameters:
% H,G,F: Channels defined as in [1].
% Mg: group size, must be small than M. Preferably a power of 2.
% varargin: algorithm parameters
%
% Output parameters:
% BDRIS: MxM BD-RIS matrix (unitary+symmetric)
% IL: Final IL value
% IL conv: convergence values.
%
% 25/03/25 -> Initial version
%
% Ignacio Santamaria, UC, 2025
%
% [1] I. Santamaria, M. Soleymani, E. Jorswieck, J. Gutierrez, "Interference 
% Minimization in Beyond-Diagonal RIS-assisted MIMO Interference Channels,"
% IEEE Open Journal of Vehicular Technology, 2025

K = size(G,2);      % number of users
M = size(G{1},2);   % number of BDRIS elements

%% Default values
opt_params = struct();
opt_params.thresholdIL = 1e-20;
opt_params.MaxiterGC = 25;
opt_params.alg = 'RtP';            % algorithm to optimize each group
if nargin < 4
    error(message('TooFewInputs'));
elseif nargin == 5
    params = varargin{1};
    for arg = fieldnames(params)'
        parameter = arg{1};
        param_value = params.(parameter);
        switch parameter
            case 'thresholdIL'
                opt_params.thresholdIL  = param_value;
            case 'MaxiterGC'
                opt_params.MaxiterGC  = param_value;
            case 'alg'
                opt_params.alg  = param_value;
        end
    end
elseif nargin > 5
    error(message('TooManyInputs'));
end

Groups = fix(M/Mg);        % number of groups
BDRIS = eye(M);            % Initialization identity matrix (this is a feasible solution)

%% Initial IL
HBDRISaux = cell(K,K);
ILconv = 0;
for i = 1:K  % tx
    for j = 1:K % rx
        HBDRISaux{i,j} = H{i,j} + F{j}*BDRIS*G{i}';  % Eq. channel for ith tx to jth rx
        if ne(i,j)
            ILconv =  ILconv + norm(HBDRISaux{i,j},'fro').^2;
        end
    end
end

%% Algorithm starts here
true = 1;
niter = 0;
while true == 1
    niter = niter+1;
    for gg = 1:Groups
        % make zero the group to update (to compute the fixed part of the
        % channel)
        BDRIS((gg-1)*Mg+1: gg*Mg,(gg-1)*Mg+1: gg*Mg) = zeros(Mg,Mg);
        % extract gth block matrix
        Fg = cell(1,K);
        Gg = cell(1,K);
        for kk=1:K
            Fg{kk} = F{kk}(:,(gg-1)*Mg+1: gg*Mg);  % block of Mg columns of F
            Gg{kk} = G{kk}(:,(gg-1)*Mg+1: gg*Mg);  % block of Mg columns of G
        end
        % Recompute the equivalent channel with the previous BDRIS block fixed
        Hg = cell(K,K);
        for i = 1:K  % tx
            for j = 1:K % rx
                Hg{i,j} = H{i,j} + F{j}*BDRIS*G{i}';  % Eq. channel for ith tx to jth rx
            end
        end
        % Compute gth block
        if strcmp(opt_params.alg,'MO')==1
            [BDRISg,~,~] = MinIL_BDRIS_MO(Hg,Gg,Fg,opt_params);
        else
            [BDRISg,~,] = MinIL_BDRIS_RtP(Hg,Gg,Fg);
        end
        BDRIS((gg-1)*Mg+1: gg*Mg,(gg-1)*Mg+1: gg*Mg) = BDRISg;
    end
    HBDRISaux = cell(K,K);
    ILg = 0;
    for i = 1:K  % tx
        for j = 1:K % rx
            HBDRISaux{i,j} = H{i,j} + F{j}*BDRIS*G{i}';  % Eq. channel for ith tx to jth rx
            if ne(i,j)
                ILg =  ILg + norm(HBDRISaux{i,j},'fro').^2;
            end
        end
    end
    ILconv = [ILconv ILg];
    %% Check convergence
    if (abs(ILconv(end)-ILconv(end-1))<opt_params.thresholdIL)||(niter >= opt_params.MaxiterGC)
        true = 0;
    end
end
IL = ILconv(end);