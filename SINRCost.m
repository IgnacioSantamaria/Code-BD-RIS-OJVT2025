function SINR = SINRCost(K,U,V,H,sigma2n)

% Description: Compute the SINR for a K-user MIMO IC with 
% fixed precoders, V, and decoders, U 
% 
% Parameters:
% K : Number of links (users)
% U : 1x K cell with decoders
% V : 1x K cell with decoders
% H : KxL cell with MIMO channels
% SINR = \sum_k || Uk^H Hkk Vk||_F^2/(\sum{k \neq l}|| Uk^H Hlk Vl||_F^2 +
% sigma2n)
%
% Ignacio Santamaria, UC 2024

Qaux = cell(K,K);
for k = 1:K
    for l = 1:K
        %------ Interference covariance matrix produced by the l-th Tx into the k-th Rx-----%
        Qaux{l,k} = U{k}'*H{l,k}*V{l}*V{l}'*H{l,k}'*U{k};
        %-----------------------------------------------------------------------------------%
    end
end
SINR = 0;         % SINR
for k = 1:K
    ind = 1:K; ind(k)=[];
    int = 0;
    for l = 1:K-1
        int = int + real(trace(Qaux{ind(l),k}));
    end
    SINR = SINR + real(trace(Qaux{k,k})/(int+sigma2n));
end

