%%------------------       Description ---------------------------------%
%  This program reproduces Fig. 6 in [1].
%  We compare the performance of group-connected BD-RIS architectures
%  optimized with MO or the RtP algorithm in [1].
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
M = [16, 24, 32, 48, 64];          % number of BD-RIS elements (Fig 6 in [1] goes to M =128, but in this case simulation takes time)
PtotaldBm = 10;                    % Transmit power in dBm (same for all users)
Ptotal = 10.^(PtotaldBm/10);       % Total power (mWatt)
B = 40;                            % Bandwidth MHz
NF = 10;                           % Noise Factor in dBs
noiseVariancedBm = -174 + 10*log10(B*10^6) + NF;
sigma2n = 10^(noiseVariancedBm/10);       % additive noise variance


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

numsimul = 10;

%% Parameters for BDRIS MO iterative algorithm
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
%INR: Interference to Noise Ratio (in dBs)
INR_noRIS = zeros(1,length(M));         % No RIS
INR_BDRISFC = zeros(1,length(M));       % iterative solution MO fully-connected
INR_RISit = zeros(1,length(M));         % RIS iterative solution (BCD)
INR_BDRIS_GC2 = zeros(1,length(M));     % Iterative relaxed Group-connected, group-size 2
INR_BDRIS_MOGC2 = zeros(1,length(M));   % Iterative MO Group-connected, group-size 2
INR_BDRIS_GC4 = zeros(1,length(M));     % Iterative relaxed Group-connected, group-size 4
INR_BDRIS_MOGC4 = zeros(1,length(M));   % Iterative MO Group-connected, group-size 4
INR_BDRIS_GC8 = zeros(1,length(M));     % Iterative relaxed Group-connected, group-size 8
INR_BDRIS_MOGC8 = zeros(1,length(M));   % Iterative MO Group-connected, group-size 8

for nn =1:numsimul   % loop for simulations (channel realizations)
    nn
    for mm = 1:length(M) % loop for BDRIS size
        mm
        %% Channel generation
        [H,G,F] = ChannelsMIMOIC(K,M(mm),Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);

        %% Group-connected Group size = 2, RtP
        opt_params.MaxiterGC = 25;         % maximum number of iterations (cycling over the groups)
        opt_params.alg = 'RtP';            % algorithm to optimize each group
        Mg = 2;
        [BDRIS2GC_RtP,IL_BDRIS2GC_RtP, IL_BDRIS2GC_RtPconv] = MinIL_GC_BDRIS(H,G,F,Mg,opt_params);

        %% Group-connected Group size = 2, MO
        opt_params.MaxiterGC = 25;         % maximum number of iterations (cycling over the groups)
        opt_params.alg = 'MO';            % algorithm to optimize each group
        Mg = 2;
        [BDRIS2GC_MO,IL_BDRIS2GC_MO, IL_BDRIS2GC_MOconv] = MinIL_GC_BDRIS(H,G,F,Mg,opt_params);

        %% Group-connected Group size = 4, RtP
        opt_params.MaxiterGC = 25;         % maximum number of iterations (cycling over the groups)
        opt_params.alg = 'RtP';            % algorithm to optimize each group
        Mg = 4;
        [BDRIS4GC_RtP,IL_BDRIS4GC_RtP, IL_BDRIS4GC_RtPconv] = MinIL_GC_BDRIS(H,G,F,Mg,opt_params);

        %% Group-connected Group size = 4, MO
        opt_params.MaxiterGC = 25;         % maximum number of iterations (cycling over the groups)
        opt_params.alg = 'MO';            % algorithm to optimize each group
        Mg = 4;
        [BDRIS4GC_MO,IL_BDRIS4GC_MO, IL_BDRIS4GC_MOconv] = MinIL_GC_BDRIS(H,G,F,Mg,opt_params);

        %% Group-connected Group size = 8, RtP
        opt_params.MaxiterGC = 25;         % maximum number of iterations (cycling over the groups)
        opt_params.alg = 'RtP';            % algorithm to optimize each group
        Mg = 8;
        [BDRIS8GC_RtP,IL_BDRIS8GC_RtP, IL_BDRIS8GC_RtPconv] = MinIL_GC_BDRIS(H,G,F,Mg,opt_params);

        %% Group-connected Group size = 8, MO
        opt_params.MaxiterGC = 25;         % maximum number of iterations (cycling over the groups)
        opt_params.alg = 'MO';            % algorithm to optimize each group
        Mg = 8;
        [BDRIS8GC_MO,IL_BDRIS8GC_MO, IL_BDRIS8GC_MOconv] = MinIL_GC_BDRIS(H,G,F,Mg,opt_params);

        %% Diagonal RIS
        opt_params.niterBCD = 1e5;            % Max number of iterations for BCD RIS opt algorithm
        [RIS,IL_RIS,ILwoRIS] = MinIL_RIS(H,G,F,opt_params);

        %% Fully-Connected BDRIS solution (MO)
        [BDRISFC_MO,IL_BDRISFC,ILBDRISFC_conv] = MinIL_BDRIS_MO(H,G,F,opt_params);

        %% INR results
        min_val = 1e-25;

        INR_BDRISFC(mm) = INR_BDRISFC(mm) + 10*log10((Ptotal*IL_BDRISFC+min_val)/sigma2n);
        INR_noRIS(mm) = INR_noRIS(mm)+ 10*log10((Ptotal*ILwoRIS+min_val)/sigma2n);
        INR_RISit(mm) = INR_RISit(mm)+10*log10((Ptotal*IL_RIS+min_val)/sigma2n);
        INR_BDRIS_GC2(mm) = INR_BDRIS_GC2(mm)+10*log10((Ptotal*IL_BDRIS2GC_RtP+min_val)/sigma2n);
        INR_BDRIS_MOGC2(mm) = INR_BDRIS_MOGC2(mm)+10*log10((Ptotal*IL_BDRIS2GC_MO+min_val)/sigma2n);
        INR_BDRIS_GC4(mm) = INR_BDRIS_GC4(mm)+10*log10((Ptotal*IL_BDRIS4GC_RtP+min_val)/sigma2n);
        INR_BDRIS_MOGC4(mm) = INR_BDRIS_MOGC4(mm)+10*log10((Ptotal*IL_BDRIS4GC_MO+min_val)/sigma2n);
        INR_BDRIS_GC8(mm) = INR_BDRIS_GC8(mm)+10*log10((Ptotal*IL_BDRIS8GC_RtP+min_val)/sigma2n);
        INR_BDRIS_MOGC8(mm) = INR_BDRIS_MOGC8(mm)+10*log10((Ptotal*IL_BDRIS8GC_MO+min_val)/sigma2n);

    end
    %plot preliminary results
    figure(2);clf; plot(M,INR_BDRISFC -INR_noRIS, 'b-pentagram','MarkerSize',ms,'LineWidth',lw);
    hold on;
    plot(M,INR_BDRIS_GC2-INR_noRIS,'b-s','MarkerSize',ms,'LineWidth',lw);
    plot(M,INR_BDRIS_MOGC2-INR_noRIS,'b--d','MarkerSize',ms,'LineWidth',lw);
    plot(M,INR_BDRIS_GC4-INR_noRIS,'r-s','MarkerSize',ms,'LineWidth',lw);
    plot(M,INR_BDRIS_MOGC4-INR_noRIS,'r--d','MarkerSize',ms,'LineWidth',lw);
    plot(M,INR_BDRIS_GC8-INR_noRIS,'g-s','MarkerSize',ms,'LineWidth',lw);
    plot(M,INR_BDRIS_MOGC8-INR_noRIS,'g--d','MarkerSize',ms,'LineWidth',lw);
    plot(M,INR_RISit-INR_noRIS,'k--*','MarkerSize',ms,'LineWidth',lw);

    legend('FC BD-RIS FC',...
        'GC BD-RIS ($M_g=2$)','GC BD-RIS MO ($M_g=2$)',...
        'GC BD-RIS ($M_g=4$)','GC BD-RIS MO ($M_g=4$)',...
        'GC BD-RIS ($M_g=8$)','GC BD-RIS MO ($M_g=8$)',...
        'RIS');
    xlabel('M (number of RIS elements)');
    ylabel('\Delta INR')
    fontname("Arial")
    fontsize(fs,"points")
    hold off;

end

INR_BDRISFC = INR_BDRISFC/numsimul;
INR_noRIS = INR_noRIS/numsimul;
INR_RISit = INR_RISit/numsimul;
INR_BDRIS_GC2 = INR_BDRIS_GC2/numsimul;
INR_BDRIS_MOGC2 = INR_BDRIS_MOGC2/numsimul;
INR_BDRIS_GC4 = INR_BDRIS_GC4/numsimul;
INR_BDRIS_MOGC4 = INR_BDRIS_MOGC4/numsimul;
INR_BDRIS_GC8 = INR_BDRIS_GC8/numsimul;
INR_BDRIS_MOGC8 = INR_BDRIS_MOGC8/numsimul;


figure(2);clf; plot(M,INR_BDRISFC -INR_noRIS, 'b-pentagram','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(M,INR_BDRIS_GC2-INR_noRIS,'b-s','MarkerSize',ms,'LineWidth',lw);
plot(M,INR_BDRIS_MOGC2-INR_noRIS,'b--d','MarkerSize',ms,'LineWidth',lw);
plot(M,INR_BDRIS_GC4-INR_noRIS,'r-s','MarkerSize',ms,'LineWidth',lw);
plot(M,INR_BDRIS_MOGC4-INR_noRIS,'r--d','MarkerSize',ms,'LineWidth',lw);
plot(M,INR_BDRIS_GC8-INR_noRIS,'g-s','MarkerSize',ms,'LineWidth',lw);
plot(M,INR_BDRIS_MOGC8-INR_noRIS,'g--d','MarkerSize',ms,'LineWidth',lw);
plot(M,INR_RISit-INR_noRIS,'k--*','MarkerSize',ms,'LineWidth',lw);

legend('FC BD-RIS FC',...
    'GC BD-RIS ($M_g=2$)','GC BD-RIS MO ($M_g=2$)',...
    'GC BD-RIS ($M_g=4$)','GC BD-RIS MO ($M_g=4$)',...
    'GC BD-RIS ($M_g=8$)','GC BD-RIS MO ($M_g=8$)',...
    'RIS');
xlabel('M (number of RIS elements)');
ylabel('\Delta INR')
fontname("Arial")
fontsize(fs,"points")
hold off;


   