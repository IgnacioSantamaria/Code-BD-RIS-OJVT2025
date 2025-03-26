function [RIS,IL,ILwoRIS] = MinIL_RIS(H,G,F,varargin)

% Description: Min-IL optimization of a BD-RIS using the Manifold
% Optimization algorithm proposed in [1]. This is Algorithm 1 in the paper
%
% Input parameters:
% H,G,F: Channels defined as in [1].
% varargin: structure with the MO parameters.
%
% Output parameters:
% BDRIS: MxM BD-RIS matrix (unitary+symmetric)
% IL: Il value after convergence
% ILBDRIS: vector with IL values vs iterations
%
% 25/03/25 -> Initial version
%
% Ignacio Santamaria, UC, 2025
%
% [1] I. Santamaria, M. Soleymani, E. Jorswieck, J. Gutierrez, "Interference 
% Minimization in Beyond-Diagonal RIS-assisted MIMO Interference Channels,"
% IEEE Open Journal of Vehicular Technology, 2025



%% Default values
opt_params = struct();
opt_params.thresholdIL = 1e-20;       % to check if IA is feasible (or to check convergence)
opt_params.niterBCD = 1e5;            % Max number of iterations for BCD RIS opt algorithm

if nargin < 3
    error(message('TooFewInputs'));
elseif nargin == 4
    params = varargin{1};
    for arg = fieldnames(params)'
        parameter = arg{1};
        param_value = params.(parameter);
        switch parameter
            case 'thresholdIL'
                opt_params.thresholdIL  = param_value;
            case 'niterBCD'
                opt_params.niterBCD  = param_value;
        end
    end
elseif nargin > 4
    error(message('TooManyInputs'));
end

niterBCD = opt_params.niterBCD;
thresholdIL = opt_params.thresholdIL;

K = size(G,2);      % number of users
M = size(G{1},2);   % number of RIS elements


%Calculate Sigma and s for RIS
Sigma = zeros(M,M);   
s = zeros(M,1);
T = 0;
for k=1:K  %rx
    for l =1:K % tx
        if ne(k,l)
            Sigma = Sigma + (G{l}.'*conj(G{l})).*(F{k}'*F{k});
            Maux = F{k}'*H{l,k}*G{l};
            s = s + diag(Maux);
            T = T + trace(H{l,k}'*H{l,k});
        end
    end
end
Sigma = (Sigma + Sigma')/2; % sanity check to make it truly psd
ILwoRIS = T;    %IL without RIS

true = 1;
iter = 0;
RIS = diag(exp(1i*2*pi*rand(M,1)));  % initial random RIS
r = diag(RIS);
ILRISpl = real(T + r'*Sigma*r + 2*real(r'*s));
while true == 1             % Alternating optimization loop
    iter = iter+1;
    theta = diag(RIS);
    for mm = 1:M            % loop to update the mth RIS element
        %mindex = 1:M;
        thetam = theta;
        %mindex(mm) = [];
        thetam(mm) = [];
        %Thetam = diag(thetam);
        sigmatilde = Sigma(:,mm); %jth column
        sigmatilde(mm) = [];
        atilde = sigmatilde'*thetam;
        theta_new = angle(s(mm) + atilde)- pi;  % new phase
        theta(mm) = exp(1i*theta_new);
    end
    RIS = diag(theta);
    r = diag(RIS);
    ILRISpl = [ILRISpl real(T + r'*Sigma*r + 2*real(r'*s))];
    
    if (abs(ILRISpl(end))<ILwoRIS*(1e-7))||(iter>=niterBCD)||((abs(ILRISpl(end)-ILRISpl(end-1))<thresholdIL))
        true = 0;
    end
end

r = diag(RIS);
% we take the abs value beacause for numerical errors we sometimes get a
% close to zero negative value
IL = abs(real(T + r'*Sigma*r + 2*real(r'*s)));

