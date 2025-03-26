function [H,G,F] = ChannelsMIMOIC(K,M,Nrx,Ntx,PosTx_XYZ,PosRx_XYZ,PosRIS_XYZ,channelparams)


% Generate K channels for the RIS-assisted
% K-user MIMO IC
%
% Input parameters:
% K : Number of users
% M: Number of RIS elements
% PosTx_XYZ, PosRx_XYZ, PosRIS_XYZ: positions, txs, rxs and RIS
% channelparams: structure with the channel parameters
%
% Output parameters
% H, G, F
%
% Ignacio Santamaria, UC 2025

H = cell(K,K);     % MIMO channels
G = cell(1,K);     % Channels from TX's to RIS
F = cell(1,K);     % Channels from RIS to RX's

pl_0 = channelparams.pl_0;
alpha_RIS = channelparams.alpha_RIS;
ray_fading = channelparams.ray_fading;
beta = channelparams.RiceRIS;
alpha_U = channelparams.alpha_direct;
beta_direct = channelparams.RiceDirect;
blocked = channelparams.blocked;

x_tx = PosTx_XYZ(:,1);
y_tx = PosTx_XYZ(:,2);
z_tx = PosTx_XYZ(:,3);
x_rx = PosRx_XYZ(:,1);
y_rx = PosRx_XYZ(:,2);
z_rx = PosRx_XYZ(:,3);
x_ris = PosRIS_XYZ(1);
y_ris = PosRIS_XYZ(2);
z_ris = PosRIS_XYZ(3);

d_RIS_rx  = sqrt((x_ris-x_rx).^2+(y_ris-y_rx).^2+(z_ris-z_rx).^2);
d_tx_RIS  = sqrt((x_ris-x_tx).^2+(y_ris-y_tx).^2+(z_ris-z_tx).^2);
pl_tx_ris_db = zeros(1,K);
pl_tx_ris_eff = zeros(1,K);
pl_ris_rx_db = zeros(1,K);
pl_ris_rx_eff = zeros(1,K);

%% Links Tx-RIS & RIS-Rx
for k = 1:K

    pl_tx_ris_db(k)  = pl_0-10*alpha_RIS*log10(d_tx_RIS(k));
    pl_tx_ris_eff(k) = 10^((pl_tx_ris_db(k))/20);

    pl_ris_rx_db(k)  = pl_0-10*alpha_RIS*log10(d_RIS_rx(k));
    pl_ris_rx_eff(k) = 10^((pl_ris_rx_db(k))/20);

    if ray_fading == 1  % if Rayleigh fading channels
        G{k} = pl_tx_ris_eff(k)* ...
            (1/sqrt(2)*randn(Ntx(k),M)+1i*1/sqrt(2)*randn(Ntx(k),M));

        F{k} = pl_ris_rx_eff(k) * ...
            (1/sqrt(2)*randn(Nrx(k),M)+1i*1/sqrt(2)*randn(Nrx(k),M));
    else
        % ======================================================
        % =============== Modeling Rician Fading ===============
        % ======================================================
        % ---- Modeling the LOS links Tx-RIS-Rx ---
        phi_AoD1=2*pi*rand;
        phi_AoA1=2*pi*rand;
        a_D_r = exp(1i*pi*(0:M-1)'*sin(phi_AoD1));
        a_D_t = exp(1i*pi*(0:Ntx(k)-1)'*sin(phi_AoA1));
        
        G{k} = pl_tx_ris_eff(k)* ...
            ((sqrt(beta)/sqrt(beta+1))*a_D_t*a_D_r' ...
            + (1/sqrt(beta+1))*(1/sqrt(2)*randn(Ntx(k),M)+1i*1/sqrt(2)*randn(Ntx(k),M)));
        % ----------------------------------------
        phi_AoD1=2*pi*rand;
        phi_AoA1=2*pi*rand;       
        a_D_r = exp(1i*pi*(0:M-1)'*sin(phi_AoD1));
        a_D_t = exp(1i*pi*(0:Nrx(k)-1)'*sin(phi_AoA1));
        F{k} = pl_ris_rx_eff(k) * ...
            ((sqrt(beta)/sqrt(beta+1))*a_D_t*a_D_r' ...
            +(1/sqrt(beta+1))*(1/sqrt(2)*randn(Nrx(k),M)+1i*1/sqrt(2)*randn(Nrx(k),M)));
    end
end


%% Direct MIMO links from ith tx to jth rx w/o RIS

d_tx_rx = zeros(K,K);
pl_tx_rx_db = zeros(K,K);
pl_tx_rx_eff = zeros(K,K);
if (ray_fading == 1)||(beta_direct==0)
    for i = 1:K
        for j = 1:K
            d_tx_rx(i,j)  = sqrt((x_tx(i)-x_rx(j))^2+(y_tx(i)-y_rx(j))^2);
            pl_tx_rx_db(i,j)  = pl_0-10*alpha_U*log10(d_tx_rx(i,j));
            pl_tx_rx_eff(i,j) = 10^((pl_tx_rx_db(i,j))/20);
            H{i,j} = pl_tx_rx_eff(i,j)*(1/sqrt(2)*randn(Nrx(j),Ntx(i))+...
                1i*1/sqrt(2)*randn(Nrx(j),Ntx(i)));
        end

    end
else  % Rician channels
    for i = 1:K
        for j = 1:K
            d_tx_rx(i,j)  = sqrt((x_tx(i)-x_rx(j))^2+(y_tx(i)-y_rx(j))^2);
            pl_tx_rx_db(i,j)  = pl_0-10*alpha_U*log10(d_tx_rx(i,j));
            pl_tx_rx_eff(i,j) = 10^((pl_tx_rx_db(i,j))/20);
            phi_AoD1=2*pi*rand;
            phi_AoA1=2*pi*rand;
            a_D_r = exp(1i*pi*(0:Nrx(j)-1)'*sin(phi_AoD1));
            a_D_t = exp(1i*pi*(0:Ntx(i)-1)'*sin(phi_AoA1));

            H{i,j} = pl_tx_rx_eff(i,j)*...
                (sqrt(beta_direct)/sqrt(beta_direct+1))*a_D_r*a_D_t' ...
                + (1/sqrt(beta_direct+1))*(1/sqrt(2)*randn(Nrx(j),Ntx(i))+1i*1/sqrt(2)*randn(Nrx(j),Ntx(i)));
        end
    end
end

if blocked == 1
    for i = 1:K
        for j = 1:K
            H{i,j} = 0* H{i,j};
        end
    end
end