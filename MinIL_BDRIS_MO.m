function [BDRIS,IL,ILBDRIS] = MinIL_BDRIS_MO(H,G,F,varargin)

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

K = size(G,2);      % number of users
M = size(G{1},2);   % number of BDRIS elements


%% Default values
opt_params = struct();
opt_params.mu = -1e6;
opt_params.niterBDRIS = 10000;
opt_params.thresholdIL = 1e-20;       %  to check convergence

if nargin < 3
    error(message('TooFewInputs'));
elseif nargin == 4
    params = varargin{1};
    for arg = fieldnames(params)'
        parameter = arg{1};
        param_value = params.(parameter);
        switch parameter
            case 'mu'
                opt_params.mu  = param_value;
            case 'niterBDRIS'
                opt_params.niterBDRIS  = param_value;
            case 'thresholdIL'
                opt_params.thresholdIL  = param_value;
        end
    end
elseif nargin > 4
    error(message('TooManyInputs'));
end

mu = opt_params.mu;
niter = opt_params.niterBDRIS;
thresholdIL = opt_params.thresholdIL;


%Calculate Sigma and s for BDRIS as in Eq. (5)
Sigma = zeros(M^2,M^2);   
s = zeros(M^2,1);
T = 0;
for k=1:K  %rx
    for l = 1:K % tx
        if ne(k,l)
            Sigma = Sigma + kron(G{l}.'*conj(G{l}),F{k}'*F{k});
            Maux = F{k}'*H{l,k}*G{l};
            s = s + Maux(:);
            T = T + trace(H{l,k}'*H{l,k});
        end
    end
end
Sigma = (Sigma + Sigma')/2; % sanity check to make it truly psd
% The IL cost function is IL= T+r'*Sigma*r +2real(r'*s)
ILwoRIS = T;

%% BDRIS passive lossless, symmetric+unitary: Iterative algorithm on the manifold 
Q = orth(randn(M,M)+1i*randn(M,M));         % Random initial point
BDRISit = Q*Q.';
A = zeros(M,M);
for k = 1:K  % rx
    for l = 1:K % tx
        if ne(k,l)
            A = A + F{k}'*H{l,k}*G{l};
        end
    end
end
r = BDRISit(:);
ILBDRIS = real(T + r'*Sigma*r + 2*real(r'*s));    % Cost function
true = 1;
iterBDRIS = 0;
while true == 1
    iterBDRIS  = iterBDRIS +1;
    Grad = zeros(M,M);
    for k = 1:K  % rx
        for l = 1:K % tx
            if ne(k,l)
                Grad = Grad + ((F{k}'*F{k})*(Q*Q.')*(G{l}'*G{l}));
            end
        end
    end
    Grad = Grad*conj(Q) + A.'*conj(Q);
    Sskew = (Q'*Grad - Grad'*Q)/2;
    Qt = Q*expm(mu*Sskew);
    Mt = Qt*Qt.';
    rt = Mt(:);
    aux = real(T + rt'*Sigma*rt + 2*real(rt'*s));
    if aux<ILBDRIS(end)  % improvement
        Q = Qt;     % update
        mu = 1.01*mu; %increase step size for next iteration
        ILBDRIS = [ILBDRIS aux];
    else               % no improvement
        mu = 0.99*mu;  % decrease step size
    end
    % Check convergence
     if length(ILBDRIS)>=2
        if (abs(ILBDRIS(end))<ILwoRIS*(1e-7))||(iterBDRIS >=niter)||(abs(mu)<1e-5)||((abs(ILBDRIS(end)-ILBDRIS(end-1))<thresholdIL))
            true = 0;
        end
    end
end
BDRIS = Q*Q.';
r = BDRIS(:);

% we take the abs value because we sometimes get a
% close to zero negative value due to numerical errors
IL = abs(real(T + r'*Sigma*r + 2*real(r'*s)));
