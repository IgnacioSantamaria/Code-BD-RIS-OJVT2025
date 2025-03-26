function cost = ILCost(K,U,V,H)

% Description: Compute the interference leakage (IL) for a K-user MIMO IC with 
% fixed precoders, V, and decoders, U 
% 
% Parameters:
% K : Number of links (users)
% U : 1x K cell with decoders
% V : 1x K cell with decoders
% H : KxL cell with MIMO channels
% cost = \sum{k \neq l}|| Uk^H Hlk Vl||_F^2
%
% Ignacio Santamaria, UC 2023

Qaux = cell(K,K);
for k = 1:K
    for l = 1:K
        %------ Interference covariance matrix produced by the l-th Tx into the k-th Rx-----%
        Qaux{l,k} = U{k}'*H{l,k}*V{l}*V{l}'*H{l,k}'*U{k};
        %-----------------------------------------------------------------------------------%
    end
end
cost = 0;         % interference leakage
for k = 1:K
    ind = 1:K; ind(k)=[];
    for l = 1:K-1
        cost = cost + real(trace(Qaux{ind(l),k}));
    end
end

