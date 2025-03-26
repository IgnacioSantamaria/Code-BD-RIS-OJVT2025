%%------------------       Description ---------------------------------%
%  This program reproduces Fig. 3 in [1].
%  We study the optimal BD-RIS deployment in the K-user MIMO-IC to
%  minimize Interference Leakage.
%  For the BDRIS optimization we use the Manifold Optimization iterative
%  algorithm (Algorithm 1 in the paper).
%
%
%  All users transmit the same power.
%
% 25/03/25 -> Initial version
%
% Ignacio Santamaria, UC, 2025
%
% [1] I. Santamaria, M. Soleymani, E. Jorswieck, J. Gutierrez, "Interference
% Minimization in Beyond-Diagonal RIS-assisted MIMO Interference Channels,"
% IEEE Open Journal of Vehicular Technology, 2025
%------------------------------------------------------------------------%

format compact;
clc
clear

%% ============== Scenario (parameters) ================
M = 64;            % number of BD-RIS elements
PtotaldBm = 10;    % Transmit power in dBm (same for all users)
Ptotal = 10.^(PtotaldBm/10); % Total power (mWatt)
B = 40;                      % Bandwidth MHz
NF = 10;                     % Noise Factor in dBs
noiseVariancedBm = -174 + 10*log10(B*10^6) + NF;
sigma2n = 10^(noiseVariancedBm/10);   % additive noise variance

%% ============== Scenario 1  ==============
K = 3;              % number of users
Ntx = [3 3 3];      % number of Tx antennas per user, must be a 1xK vector
Nrx = [3 3 3];      % number of Rx antennas per user, must be a 1xK vector

%% ============== Channel parameters ==============
channelparams = struct;
channelparams.blocked = 0;         % Set to 1 if direct channel is blocked
channelparams.RiceRIS = 3;         % Rician factor for the Tx-RIS-Rx channels 
channelparams.RiceDirect = 0;      % Rician factor for the direct links (if 0 the fading is Rayleigh)
channelparams.pl_0 = -28;          % Path loss at a reference distance (d_0)
channelparams.alpha_RIS = 2;       % Path loss exponent for the RIS links
channelparams.alpha_direct = 3.75; % Path loss exponent for the direct links
channelparams.ray_fading = 0;      % Set to 1 if all channels Rayleigh

%% --- Position of the users/RIS (units in meters)-----%%
sqr_size = 50;
indices = 0:1/(K-1):1;
% At the x-axis
x_tx = zeros(1,K);
x_rx = sqr_size + x_tx;
x_ris = 5:5:45;  % we vary the x coordinate for the BDRIS in 5 meters steps

% At the y-axis
y_tx=sqr_size*indices;   % users esquispaced along the y-axis
y_rx = y_tx;
y_ris = 5:5:45;   % we vary the y coordinate for the BDRIS in 5 meters steps

% At the z-axis
z_tx  = 1.5*ones(1,K);
z_rx  = z_tx;
z_ris = 5;


PosTx_XYZ = [x_tx' y_tx' z_tx'];
PosRx_XYZ = [x_rx' y_rx' z_rx'];

numsimul = 2;  % For improved results, increase the number of simulations

%% Parameters for BDRIS/RIS iterative algorithms
opt_params.mu = -1e6;
opt_params.niterBDRIS = 2000;
opt_params.thresholdIL = 1e-20;       % to check convergence

%% Parameters for figures
fs = 12;   % fontsize
lw = 1.5;  % linewidth
ms = 8;    % markersize
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Store variables
%INR: Interference to Noise Ratio (in dBs) (it is worst-case INR for any unitary precoder/decoder)
INR_noRIS = zeros(length(y_ris),length(x_ris));     % No RIS
INR_BDRIS = zeros(length(y_ris),length(x_ris));     % MO solution (iterative)
HBDRIS = cell(K,K);

for nn =1:numsimul
    nn
    for xx = 1:length(x_ris)
        for yy = 1:length(y_ris)
            disp(['Position x = ' num2str(xx), ',' 'Position y = ' num2str(yy)])
            %% Channel generation
            PosRIS_XYZ = [x_ris(xx)', y_ris(yy)', z_ris'];
            [H,G,F] = ChannelsMIMOIC(K,M,Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);
            
            %% BDRIS solution (iterative algorithm)
            %disp('Start BDRIS iterative algorithm')
            [BDRIS,ILpl_it,ILBDRISconv] = MinIL_BDRIS_MO(H,G,F,opt_params);
            %disp('Finished passive lossless BDRIS')
            INR_noRISaux = 0;
            for i = 1:K  % tx
                for j = 1:K % rx
                    HBDRIS{i,j} = H{i,j} + F{j}*BDRIS*G{i}';  % Equivalent channel for ith tx to jth rx
                    if ne(i,j)
                        INR_noRISaux =  INR_noRISaux + norm(H{i,j},'fro').^2;
                    end
                end
            end
            figure(1);plot(ILBDRISconv);
            IL_BDRIS = ILpl_it;
            %%INR results
            min_val = 1e-25;  % to avoid numerical issues
            INR_BDRIS(yy,xx) = INR_BDRIS(yy,xx) + 10*log10((Ptotal*IL_BDRIS+min_val)/sigma2n);
            INR_noRIS(yy,xx) = INR_noRIS(yy,xx) + 10*log10((Ptotal*INR_noRISaux+min_val)/sigma2n);
        end
    end
end
INR_BDRIS = INR_BDRIS/numsimul;
INR_noRIS = INR_noRIS/numsimul;
figure(1);contour(y_ris,x_ris,INR_BDRIS - mean(mean(INR_noRIS)),6,'ShowText','on',"LabelFormat","%0.1f")
figure(2);contourf(y_ris,x_ris,INR_BDRIS-mean(mean(INR_noRIS)),6,'ShowText','on',"LabelFormat","%0.1f")


% The optimal coordinates for Scenario 1 are roughly (40,25) (x[m], y[m])