function [U,V,BDRIS,ILouter,ILinner] = MinIL_BDRISJointDesign(H,G,F,U,V,varargin)

% Desc: Obtains precoders, decoders and a BDRIS (symmetric,
% globally passive) that minimize the interference leakage
%
% Parameters
% H : KxK cell with MIMO channels
% G : 1xK cell with channels from Tx's to RIS
% F : 1xK cell with channels from RIS to RX's
% d : 1xK vector with number of streams per user
% niterAO:  maximum number iteration alternating optimization (if needed)
% thresholdIL: Small IL value to check convergence


% U : 1x K cell with decoders
% V : 1xL cell with precoders
% RIS : MxM diagonal RIS matrix

K = size(G,2);      % number of users
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
opt_params.thresholdinner2 = 1e-3;  % To check convergence of the inner loop (decrease between two consecutive iterations)
opt_params.niterouter = 200;        % Maximum number of iterations for outer loop (RIS optimization)
opt_params.thresholdouter1 = 1e-8;  % To check convergence of the outer loop (absolute ISP value)
opt_params.thresholdouter2 = 1e-3;  % To check convergence of the outer loop (decrease between two consecutive iterations)

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
            case 'niterouter'
                opt_params.niterouter  = param_value;
            case 'thresholdouter1'
                opt_params.thresholdouter1  = param_value;
            case 'thresholdouter2'
                opt_params.thresholdouter2  = param_value;
        end
    end
elseif nargin > 6
    error(message('TooManyInputs'));
end

ILwoRIS = ILCost(K,U,V,H);    % initial IL without RIS

true_outer = 1;
nn_outer = 0;
ILouter = zeros(1,opt_params.niterouter);    % Stores interference leakage
%% Outer loop starts here
while (true_outer==1)

    nn_outer = nn_outer+1;
    %% -------------------------------  BDRIS update --------------------------------%
    [BDRIS,~] = MinIL_BDRISsymgp(H,G,F,U,V);  % optimize BDRIS
    %% Equivalent channel
    Heq = cell(K,K);
    for i = 1:K  % tx
        for j = 1:K % rx
            Heq{i,j} = H{i,j} + F{j}*BDRIS*G{i}';  % Heq(i,j) channel for ith tx to jth rx with RIS
        end
    end
    ILouter(nn_outer) = ILCost(K,U,V,Heq);
    %% Check outer convergence outer loop
    if nn_outer==1
        if ILouter(nn_outer)<ILwoRIS*opt_params.thresholdouter1
            true_outer = 0;
            ILouter(nn_outer+1:end) = [];
            ILinner = ILouter;
        end
    else
        if (abs(ILouter(nn_outer)-ILouter(nn_outer-1))<ILouter(nn_outer-1)*opt_params.thresholdouter2)||...
                (ILouter(nn_outer)<ILwoRIS*opt_params.thresholdouter1)...
                ||(nn_outer==opt_params.niterouter)
            true_outer = 0;
            ILouter(nn_outer+1:end) = [];
        end
    end
    if true_outer == 1
        [U,V,ILinner] = MinIL_UV(Heq,ILwoRIS,U,V,opt_params);
    end

end


