function [U,V,SINR] = MaxSINR_UV(Heq,U,V,sigma2n,Ptotal,varargin)

% Precoder and decoder optimization for SINR maximization
% for the K-user MIMO-IC. Stage II

K = size(V,2);      % number of users
Ntx = zeros(1,K);   % transmit antennas
Nrx = zeros(1,K);   % receive antennas
d = zeros(1,K);     % streams
for kk = 1:K
    Ntx(kk) = size(V{kk},1);
    Nrx(kk) = size(U{kk},1);
    d(kk) = size(V{kk},2);
end
%% Default values
opt_params = struct();
opt_params.niterinner = 100;        % Maximum number of iterations for inner loop (U&V opt)
opt_params.thresholdinner1 = 1e-8;  % To check convergence of the inner loop (absolute  value)
opt_params.thresholdinner2 = 1e-4;  % To check convergence of the inner loop (increase between two consecutive iterations)
if nargin < 5
    error(message('TooFewInputs'));
elseif nargin == 6
    params = varargin{1};
    for arg = fieldnames(params)'
        parameter = arg{1};
        param_value = params.(parameter);
        switch parameter
            case 'niterinner'
                opt_params.niterinner  = param_value;
            case 'thresholdinner1'
                opt_params.thresholdinner1  = param_value;
            case 'thresholdinner2'
                opt_params.thresholdinner2  = param_value;
        end
    end
elseif nargin > 6
    error(message('TooManyInputs'));
end

%% -------------------------------  U&V update --------------------------------%
true_inner = 1;
nn_inner = 0;
SINR = zeros(1,opt_params.niterinner);    % Stores SINR
%SINRini = SINRCost(K,U,V,H,sigma2n);
while true_inner    % here Heq does not change
    nn_inner = nn_inner+1;

    %% -------------------------- Decoder update (U) ---------------------------------%
    Qaux = cell(K,K);       % interference covariance matrices (K^2 matrices)
    Q = cell(1,K);          % interference covariance matrices at the k-th Rx (K matrices)
    for k = 1:K
        for l = 1:K
            %-- Interference covariance matrix produced by the l-th Tx into the k-th Rx----%
            Qaux{l,k} = Heq{l,k}*V{l}*V{l}'*Heq{l,k}';
            %------------------------------------------------------------------------------%
        end
    end
    for k = 1:K
        ind = 1:K; ind(k)=[];   % a vector with all indexes except the k-th
        %---- Int. convariance matrix at the k-th Rx ---%
        Q{k} = zeros(size(Qaux{k,k}));
        for tt =1:length(ind)
            Q{k} = Q{k} + Qaux{ind(tt),k};
        end
        %------ Now we compute the the interference subspace and  ----%
        %----- its orthogonal complement------------------------------%
        %DirectChannel = Heq{k,k}*V{k}*V{k}'*Heq{k,k}';
        DirectChannel = Heq{k,k}*V{k};
        %[A,~] = svd((Q{k}+sigma2n*eye(size(DirectChannel,1)))\DirectChannel);
        A = orth((Q{k}+sigma2n*eye(size(DirectChannel,1)))\DirectChannel);
        U{k} = A(:,1:d(k));  % largest eigenvectors -> interference free subspace-> Projector U(:,:,k)*U(:,:,k)'
    end


    %% -------------------------------  Precoders update (V) --------------------------------%
    Qaux = cell(K,K);    % interference covariance matrices (K^2 matrices)
    Q = cell(1,K);       % interference covariance matrices provoked by the l-th Tx (K matrices)
    for k = 1:K
        for l = 1:K
            %------ Interference covariance matrix seen by the l-th Tx from the k-th Rx-----%
            Qaux{l,k} = Heq{l,k}'*U{k}*U{k}'*Heq{l,k};
            %--------------------------------------------------------------
        end
    end

    for l = 1:K
        ind = 1:K;
        ind(l) = [];       % we take out the desired user
        Q{l} = zeros(size(Qaux{l,l}));
        for tt = 1:length(ind)
            Q{l} = Q{l} + Qaux{l,ind(tt)};
        end
        %------ find the subspace that maximizes SINR ----%
        %DirectChannel = Heq{l,l}'*U{l}*U{l}'*Heq{l,l};
        DirectChannel = Heq{l,l}'*U{l};
        %[A,~] = svd((Q{l}+sigma2n*eye(size(DirectChannel,1)))\DirectChannel);
        A = orth((Q{l}+sigma2n*eye(size(DirectChannel,1)))\DirectChannel);
        V{l} = A(:,1:d(l)); % smallest eigenvectors of the interference subspace
    end
    for kk = 1:K
            V{kk} = sqrt(Ptotal)*(V{kk}/sqrt(trace(V{kk}'*V{kk})));
    end
    %% IL Cost
    SINR(nn_inner) = SINRCost(K,U,V,Heq,sigma2n);
    %% Check convergence
    if nn_inner>1
        
        if ((SINR(nn_inner)-SINR(nn_inner-1))<SINR(nn_inner-1)*opt_params.thresholdinner2)||...
               (nn_inner==opt_params.niterinner)
            true_inner = 0;
            SINR(nn_inner+1:end) = [];
        end
    end
end