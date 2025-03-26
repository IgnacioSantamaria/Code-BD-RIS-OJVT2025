function C = function_sumcap(H,V,var,K)

%%----------------------- Description ------------------------------------%
% Obtains the user rates for an interference MIMO channel
% The rate for each user is calculated as
% C(k)=log2(abs(det((eye(N)+Ik)\Qk)), where Ik is the NxN total interference
% matrix and Qk is the transmitted covariance matrix for the k-th user user
%
% H: Kx K cell with (equivalent) MIMO Int. channels
% V: 1x K cell with precoders 
% var: Noise variance
% K: number of users
% C: Kx1 vector with the user rates
%
% Ignacio Santamaria, UT Austin, Dec. 09
% Rev. Feb. 2024
% Modified Sept. 24

C = zeros(K,1);          % store the achievable rates for each user
Qaux = cell(K,K);        % interference covariance matrices (K^2 matrices)
Q = cell(1,K);            % interference covariance matrices at the k-th Rx (K matrices)
for k = 1:K
    for l = 1:K
        %-- Interference covariance matrix produced by the l-th Tx into the k-th Rx----%
        Qaux{l,k} = H{l,k}*V{l}*V{l}'*H{l,k}';
        %------------------------------------------------------------------------------%
    end
end
for k = 1:K
    ind = 1:K; ind(k) = [];   % a vector with all indexes except the k-th
    Q{k} = zeros(size(Qaux{k,k}));
    for tt = 1:length(ind)
        Q{k} = Q{k} + Qaux{ind(tt),k};
    end
    %---- Int. + noise covariance matrix at the k-th Rx ---%
    Q{k} = Q{k} + var*eye(size(Q{k},2));
    C(k) = log2(abs(det(eye(size(Q{k},2))+ Q{k}\Qaux{k,k})));
end

