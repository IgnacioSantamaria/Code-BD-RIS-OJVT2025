%%------------------       Description ---------------------------------%
%  This program compares the manifold optimization algorithm vs the RtP
%  algorithm. Algorithms 1 & 2 in [2].
%  For small values of M (M<16) the RtP algorithm provides a reasonably good approximation
%  For M>16 the RtP solution might be very suboptimal. 
%
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
M = 8;             % number of BD-RIS elements
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
x_ris = 40;  % Optimal x coordinate for BD-RIS

% At the y-axis
y_tx=sqr_size*indices;   % users esquispaced along the y-axis
y_rx = y_tx;
y_ris = 25;   % optimal t coordinate for BD-RIS

% At the z-axis
z_tx  = 1.5*ones(1,K);
z_rx  = z_tx;
z_ris = 5;

PosTx_XYZ = [x_tx' y_tx' z_tx'];
PosRx_XYZ = [x_rx' y_rx' z_rx'];
PosRIS_XYZ = [x_ris', y_ris', z_ris'];

numsimul = 1;

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

%% Channel generation
[H,G,F] = ChannelsMIMOIC(K,M,Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);
%% IL without BD-RIS
ILwoRIS = 0;
for i = 1:K  % tx
    for j = 1:K % rx
        if ne(i,j)
            ILwoRIS =  ILwoRIS + norm(H{i,j},'fro').^2;
        end
    end
end

for nn = 1:numsimul   % loop for simulations (MO runs starting from different initialization points)
    nn
    %% BDRIS solution (MO iterative algorithm)
    [BDRIS_MO,IL_MO,ILBDRISconvMO] = MinIL_BDRIS_MO(H,G,F,opt_params);
  
    %% BDRIS solution (RtP solution)
    [BDRIS_RtP,IL_RtP] = MinIL_BDRIS_RtP(H,G,F);

    %% Polt results
    a = 10*log10(ILBDRISconvMO/ILwoRIS);
    figure(1);clf;
    plot(a,'r');
    hold on;
    plot(10*log10(IL_RtP/ILwoRIS)*ones(size(a)),'b');
    xlabel('Manifold Optimization (MO) Iterations')
    ylabel('\Delta INR (dB)')
    legend('MO', 'RtP');
    hold off
end





