function [U,V,IL] = MinIL_UV(Heq,ILwoRIS,U,V,varargin)

% Precoder and decoder optimization for IL minimization
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
opt_params.thresholdinner1 = 1e-8;  % To check convergence of the inner loop (absolute ISP value)
opt_params.thresholdinner2 = 1e-4;  % To check convergence of the inner loop (decrease between two consecutive iterations)
if nargin < 4
    error(message('TooFewInputs'));
elseif nargin == 5
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
elseif nargin > 5
    error(message('TooManyInputs'));
end

%% -------------------------------  U&V update --------------------------------%
true_inner = 1;
nn_inner = 0;
IL = zeros(1,opt_params.niterinner);    % Stores interference leakage
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
        [A,~] = svd(Q{k});
        U{k} = A(:,Nrx(k)-d(k)+1:end);  % smallest eigenvectors -> interference free subspace-> Projector U(:,:,k)*U(:,:,k)'
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
        %------ find the subspace that causes less interference ----%
        [A,~] = svd(Q{l});
        V{l} = A(:,Ntx(l)-d(l)+1:end); % smallest eigenvectors of the interference subspace
    end
    %% IL Cost
    IL(nn_inner) = ILCost(K,U,V,Heq);
    %% Check convergence
    if nn_inner==1
        if IL(nn_inner)<ILwoRIS*opt_params.thresholdinner1
            true_inner = 0;
            IL(nn_inner+1:end) = [];
        end
    else
        if (abs(IL(nn_inner)-IL(nn_inner-1))<IL(nn_inner-1)*opt_params.thresholdinner2)||...
                (IL(nn_inner)<ILwoRIS*opt_params.thresholdinner1)...
                ||(nn_inner==opt_params.niterinner)
            true_inner = 0;
            IL(nn_inner+1:end) = [];
        end
    end
end