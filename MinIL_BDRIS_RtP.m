function [BDRIS,IL] = MinIL_BDRIS_RtP(H,G,F)

% Description: Min-IL optimization of a BD-RIS using the Relax-then-Project 
% (RtP) algorithm proposed in [1]. This is Algorithm 2 in the paper
%
% Input parameters:
% H,G,F: Channels defined as in [1].
%
% Output parameters:
% BDRIS: MxM BD-RIS matrix (unitary+symmetric)
% IL: Final IL value 
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

%% Calculate T, Sigma and s for a BD-RIS (as in Eq. (5) in the paper)
Sigma = zeros(M^2,M^2);
s = zeros(M^2,1);
T = 0;
for k=1:K  %rx
    for l =1:K % tx
        if ne(k,l)
            Sigma = Sigma + kron(G{l}.'*conj(G{l}),F{k}'*F{k});
            Maux = F{k}'*H{l,k}*G{l};
            s = s + Maux(:);
            T = T + trace(H{l,k}'*H{l,k});
        end
    end
end
Sigma = (Sigma + Sigma')/2; % sanity check to make it truly psd (not really needed)
% The IL cost function is IL = T + r'*Sigma*r + 2real(r'*s)

A = reshape(1:M*M, M, M);
v = reshape(A', 1, []);
P = eye(M^2);
P = P(v,:);    % This is the conmutation matrix
               % used for transforming the vectorized form of a matrix into the vectorized form of its transpose

SymmetryLC =(P-eye(M^2));  % The solution must belong to the nullity of this matrix,
rankSymmetryLC = M*(M-1)/2;
[Uslc,~,~] = svd(SymmetryLC);
NoiseSubs = Uslc(:,rankSymmetryLC+1:end);   % The solution must belong to the nullity of the subspace

Sigmasymm = NoiseSubs'*Sigma*NoiseSubs;             % here we rename Sigma!
ssymm = NoiseSubs'*s;                               % here we rename s!

% Problem reformulation in terms of r (Eq. (12) in [1])
[USigma,D,~] = svd(Sigmasymm);  % This is still a comp. bottleneck!!
d = diag(D);
SigmaRank = length(find(d>max(size(D))*eps(norm(D))));  % This is the rank of Sigma
Ds = D(1:SigmaRank,1:SigmaRank);
ds = diag(Ds);
Us = USigma(:,1:SigmaRank);        % Basis for the signal subspace
sprima = Us'*ssymm;
rs_star = -diag(1./ds)*sprima;
r_unc = Us*rs_star;
norm_unc = real(r_unc'*r_unc);
if norm_unc>M  % here we use bisection
    %disp('norm_unc >M')
    % we need to find a lambda s.t the norm of the solution is smaller than M
    lambda_min = 0;
    norm_max = norm_unc;
    lambda_max = 1e-11;
    true = 1;
    while true==1
        lambda_max =  lambda_max*10;
        rs_star = -diag(1./(ds+lambda_max))*sprima;
        r_star = Us*rs_star;
        norm_min = real(r_star'*r_star);
        if norm_min<M
            true = 0;
        end
    end
    %% Start bisection
    true_bis = 1;
    while true_bis == 1
        lambda = (lambda_min + lambda_max)/2;
        rs_star = -diag(1./(ds+lambda))*sprima;
        r_star = Us*rs_star;
        norm_lambda = real(r_star'*r_star);
        if norm_lambda<M
            lambda_max = lambda;
            norm_min = norm_lambda;
        elseif norm_lambda>M
            lambda_min = lambda;
            norm_max = norm_lambda;
        elseif norm_lambda==M
            lambda_min = lambda;
            lambda_max = lambda;
            norm_min = M;
            norm_max = M;
        end
        %Check covergence 
        if (abs(norm_min-norm_max)<M*1e-5)
            true_bis = 0;
        end
    end
    % The final solution is r_star

else  
    %disp('norm_unc<M')
    r_star = r_unc;
end


r_star = NoiseSubs*r_star; % This is now an M^2x1 vector (see Eq 12)
BDRISrelax = reshape(r_star,M,M);

%% Projection onto the set of unitary matrices (Takagi)

BDRIS_aux = (BDRISrelax + BDRISrelax.')/2;  %This is not needed (already symmetric) just a sanity check
[Uaux,Daux,Vaux] = svd(BDRIS_aux);
daux= diag(Daux);
rangobis = length(find(daux>max(size(Daux))*eps(norm(Daux))));
Uaux = [Uaux(:,1:rangobis) conj(Vaux(:,rangobis+1:end))];
BDRIS = Uaux*Vaux';          % unitary + symmetric
r_star = BDRIS(:);           % Final passive+symmetric solution
IL = abs(real(T + r_star'*Sigma*r_star + 2*real(r_star'*s)));    %Final IL
