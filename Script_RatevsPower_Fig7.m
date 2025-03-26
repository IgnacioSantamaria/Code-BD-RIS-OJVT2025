%%------------------       Description ---------------------------------%
%  This program reproduces Fig. 7 in [1].
%  We compare the rate vs transmit power performance of BD-RIS-assisted
%  MIMO ICs
%  
% In Stage I the BD-RIS is optimized to minimize IL using the MO algorithm
% In Stage II the precoders (and decoders) are optimized according to
% different criteria, namely: SVD-based precoders, min-IL precoders, and
% max-SINR precoders. These are a subset of the precoders compread in Fig.
% 6 of [1]
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
M = 16;                            % number of BD-RIS elements 
PtotaldBm = -4:4:32;               % Transmit power in dBm (same for all users)
Ptotal = 10.^(PtotaldBm/10);       % Total power (mWatt)
B = 40;                            % Bandwidth MHz
NF = 10;                           % Noise Factor in dBs
noiseVariancedBm = -174 + 10*log10(B*10^6) + NF;
sigma2n = 10^(noiseVariancedBm/10);       % additive noise variance

%% ============== Scenario 1  ==============
K = 3;              % number of users
Ntx = [3 3 3];      % number of Tx antennas per user, must be a 1xK vector
Nrx = [3 3 3];      % number of Rx antennas per user, must be a 1xK vector
d = [2 2 2];        % number of streams transmitted per user, must be a 1xK vector d(i) <= min(Ntx(i),Nrx(i))

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

%% Num Simulations
numsimul = 25;      % Increase the number of simulations for smoother curves.

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
SR_noRIS = zeros(1,length(PtotaldBm));         % No RIS
SR_SVD = zeros(1,length(PtotaldBm));           % SVD precoders (disregarding residual interference)
SR_minIL = zeros(1,length(PtotaldBm));         % Min-IL precoders 
SR_maxSINR = zeros(1,length(PtotaldBm));       % Max-SINR precoders 

%Equivalent channels
HBDRISeq = cell(K,K);
HBDRISjoint = cell(K,K);

for nn =1:numsimul   % loop for simulations (channel realizations)
    nn
    for mm = 1:length(PtotaldBm) % loop for BDRIS size
        mm
        P=10.^(PtotaldBm(mm)/10);
        %% Channel generation
        [H,G,F] = ChannelsMIMOIC(K,M,Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);

        ILwoRIS = 0;  % Int. Leakage without BDRIS
        for i = 1:K  % tx
            for j = 1:K % rx
                if ne(i,j)
                    ILwoRIS = ILwoRIS+ norm(H{i,j},'fro')^2; % IL
                end
            end
        end       

        %% 1st stage min-IL BD-RIS optimized with MO algorithm
        if ne(M,0)
            %% BDRIS solution (MO iterative algorithm)
            [BDRIS_MO,IL_MO,ILBDRISconvMO] = MinIL_BDRIS_MO(H,G,F,opt_params);
        else
            BDRIS_MO = 0;  % No BDRIS
        end

        %% 2nd stage precoder/decoder design

        %===== SVD-based precoder/decoder (neglecting interference, non-cooperative design)======%
        % Find equivalent channel
        for i = 1:K  % tx
            for j = 1:K % rx
                HBDRISeq{i,j} = H{i,j} + F{j}*BDRIS_MO*G{i}';  % Eq. channel for ith tx to jth rx
            end
        end
        % Precoders & decoders
        Vsvd = cell(1,K);   % Precoders
        Usvd = cell(1,K);   % Decoders

        for kk = 1:K
            [Qaux,Vsvd{kk},paux,Usvd{kk}] = OptTransmitCovMatrix(HBDRISeq{kk,kk},sigma2n*eye(Nrx(kk)),Ptotal(mm));
            indp = find(paux~=0);  % eliminate zeros. Indices with zero allocated power
            paux = paux(indp);
            Vsvd{kk} = Vsvd{kk}(:,indp);
            Vsvd{kk} = Vsvd{kk}*diag(sqrt(paux));  % with power allocation
        end

        %% ==================== Min-IL design ===================%
        % Generate initial precoders & decoders
        VminIL = cell(1,K);   % Precoders
        UminIL = cell(1,K);   % Decoders
        for kk = 1:K
            VminIL{kk} = sqrt(Ptotal(mm)/Ntx(kk))*orth(randn(Ntx(kk),d(kk))+1i*randn(Ntx(kk),d(kk)));
            UminIL{kk} = orth(randn(Nrx(kk),d(kk))+1i*randn(Nrx(kk),d(kk)));
        end
        [UminIL,VminIL,IL_after_MinIL_prec] = MinIL_UV(HBDRISeq,ILwoRIS,UminIL,VminIL);
        %Power normalization
        for kk=1:K
            VminIL{kk} = sqrt(Ptotal(mm))*(VminIL{kk}/sqrt(trace(VminIL{kk}'*VminIL{kk})));
        end

        %% Max SINR design
        VmaxSINR = cell(K,K);
        UmaxSINR = cell(K,K);
        for kk = 1:K
            VmaxSINR{kk} = sqrt(Ptotal(mm)/Ntx(kk))*orth(randn(Ntx(kk),d(kk))+1i*randn(Ntx(kk),d(kk)));
            UmaxSINR{kk} = orth(randn(Nrx(kk),d(kk))+1i*randn(Nrx(kk),d(kk)));
        end
        [UmaxSINR,VmaxSINR,SINR] = MaxSINR_UV(HBDRISeq,UmaxSINR,VmaxSINR,sigma2n,Ptotal(mm),opt_params);
        %Power normalization
        for kk=1:K
            VmaxSINR{kk} = sqrt(Ptotal(mm))*(VmaxSINR{kk}/sqrt(trace(VmaxSINR{kk}'*VmaxSINR{kk})));
        end

        
        %% Results
        %Sum rate results
        SR_SVD(mm) = SR_SVD(mm) + mean(function_sumcap(HBDRISeq,Vsvd,sigma2n,K));
        SR_minIL(mm) = SR_minIL(mm) + mean(function_sumcap(HBDRISeq,VminIL,sigma2n,K));
        SR_maxSINR(mm) = SR_maxSINR(mm) + mean(function_sumcap(HBDRISeq,VmaxSINR,sigma2n,K));  
        

    end
    figure(1);plot(PtotaldBm,K*SR_SVD, 'b-pentagram','MarkerSize',ms,'LineWidth',lw);
    hold on;
    plot(PtotaldBm,K*SR_minIL, 'r-o','MarkerSize',ms,'LineWidth',lw);
    plot(PtotaldBm,K*SR_maxSINR, 'g-d','MarkerSize',ms,'LineWidth',lw);
    legend('SVD', 'Min-IL', 'Max-SINR');
    xlabel('P_t (dB)');
    ylabel('Sum rate (b/s/Hz)')
    fontname("Arial")
    fontsize(fs,"points")
    hold off

end
SR_SVD = SR_SVD/numsimul;
SR_minIL = SR_minIL/numsimul;
SR_maxSINR = SR_maxSINR/numsimul;

figure(20);clf; plot(PtotaldBm,K*SR_SVD, 'b-pentagram','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(PtotaldBm,K*SR_minIL, 'r-o','MarkerSize',ms,'LineWidth',lw);
plot(PtotaldBm,K*SR_maxSINR, 'g-d','MarkerSize',ms,'LineWidth',lw);
legend('SVD', 'Min-IL', 'Max-SINR');
xlabel('P_t (dB)');
ylabel('Sum rate (b/s/Hz)')
fontname("Arial")
fontsize(fs,"points")
hold off;