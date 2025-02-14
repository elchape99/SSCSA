clc
close all
clear 
% Script developed for use two piezo

csi_i = 0; % natural damping of the system
% with different damping change the peak, but the intersection point with
% the short circuit response does not change
csi1 = 0.005;  
csi2 = 0.005;  
k = [0.24, 0.18];
k1= k(1);
k2= k(2);

% 
%w_i_oc_oc = [];
% w_i: freq naturali

% w_1 = w_i(1);
% w_2 = w_i(2);
w_1 = 0;%20.5;
w_2 = 116;
% w_1: freq ottime piezo 1
% w_2: freq ottime piezo 2

% data of the bar 
beam_defm;

% w   = (1:0.5:500)';
% w_i = [20.122, 115.74]; % w_sc_sc
% w_cap = [22.723,117.7]; % w_oc_oc

w  = 2*pi* (1:0.1:500)';
w_i = 2*pi* beam.nf.sc_sc(1:2); % w_sc_sc
w_cap = 2*pi* beam.nf.oc_oc(1:2); % w_oc_oc

%optimal values if uncorrelated
w_e_1_opt = w_cap(1);
% w_e_1_opt = w_i(1);
% w_e_1_opt = 0;
w_e_2_opt = w_cap(2);

csi_e_1_opt = sqrt(3)/2 * sqrt ((w_cap(1)^2 - w_i(1)^2) / (w_cap(1)^2 + w_i(1)^2));
csi_e_2_opt = sqrt(3)/2 * sqrt ((w_cap(2)^2 - w_i(2)^2) / (w_cap(2)^2 + w_i(2)^2));

L1_opt = 1/ (w_e_1_opt^2 * beam.Cp.C11);
L2_opt = 1/ (w_e_2_opt^2 * beam.Cp.C22);

R1_opt = 2 * csi_e_1_opt / (beam.Cp.C11 * w_e_1_opt);
R2_opt = 2 * csi_e_2_opt / (beam.Cp.C22 * w_e_2_opt);

% L1_opt = inf;
% R1_opt = inf;


% H_rl_rl = double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], [beam.Cp.C21, beam.Cp.C22], L1, L2, R1, R2, beam.k.k1(1:2), beam.k.k2(1:2));
H_rl_rl = double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], [beam.Cp.C21, beam.Cp.C22], L1_opt, L2_opt, R1_opt, R2_opt, beam.k.k1(1:2), beam.k.k2(1:2));


% H_rl_rl = ((-w.^2 + 2 .* 1i .* csi1 .* w_1 .* w + w_1.^2) .* (-w.^2 + 2 .* 1i .* csi2 .* w_2 .* w + w_2.^2)) ...
%     ./ ...
%     (-w.^6 + ...
%     2.* 1i .* w.^5 .* (csi1 .* w_1 + csi2 .* w_2 + csi_i .* w_i) ...
%     + w.^4 .* (w_2.^2 + w_1.^2 + 4 .* csi1 .* csi2 .* w_1 .* w_2 + 4 .* csi_i .* w_i .* (csi1 .* w_1 + csi2 .* w_2) + w_cap.^2) ...
%     - 2 .* 1i .* w.^3 .* (csi1 .* w_1 .* w_2.^2 + csi2 .* w_2 .* w_1.^2 + csi_i .* w_i .* (w_2.^2 + w_1.^2 + 4 .* csi1 .* csi2 .* w_1 .* w_2) + w_cap.^2 .* (csi1.* w_1 + csi2 .* w_2))...
%     - w.^2 .* (w_1.^2 .* w_2.^2 + 4 .* csi_i .* w_i .* (csi1 .* w_1 .* w_2.^2 + csi2 .* w_2 .* w_1.^2) + w_cap.^2 .* (w_1.^2 + w_2.^2 + 4 .* csi1.*csi2 .* w_1 .* w_2) ...
%     + 2 .* 1i .* w .* (csi_i .* w_i .* w_1.^2 .* w_2.^2 + w_cap.^2 .* (csi1 .* w_1 .* w_2.^2 + csi2 .* w_2 .* w_1^2)) ...
%     + w_cap.^2 .* (w_1.^2 .* w_2.^2)));

j = 1i;

% H_rl_rl_num = ((-w.^2 + 2*j.*csi1.*w_1 + w_1.^2).*(-w.^2 + 2*j.*csi2.*w_2 + w_2.^2));
        
% H_rl_rl_den = (-w.^6 + ...
%                 2jw^5 .* (csi1w_1 + csi2w_2 + csi_iw_i) + ...
%                 w^4 .* (w_1^2 + w_2^2 + 4csi1csi2w_1w_2 + 4csi_iw_i.*(csi1w_1 + csi2w_2) + w_cap^2) ...
%                 -2jw^3)

      
      
H_sc_sc = 1./(w_i.^2 - w.^2 + 1i.*2.*w_i.*w.*csi_i);

% plotting the frf
figure(1)
% semilogy(w, abs(sum(H_rl_rl, 2)), w, abs(sum(H_sc_sc,2)));
% legend("H-rl-rl", "H-sc-sc")
semilogy(w, abs(H_rl_rl), w, abs(sum(H_sc_sc,2)));
legend("H-rl-rl-1", "H-rl-rl-2", "H-sc-sc")

FRF_rl_rl_double_piezo_first = load ("FRF_rl_rl_double_piezo_first");
FRF.sc_sc = load ("FRF_sc_sc.mat");

% x comumn for the plot
x = linspace(1,500, length(FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl));

% plot the total FRF, all processing for mode 1 (sc, oc, only pezo 1)
figure
semilogy (x, FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl, "LineWidth",0.6);
hold on
semilogy(x, FRF_rl_rl_double_piezo_first.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl, "LineWidth",0.6);
hold on
title ("FRF-sc-sc - FRF-oc-sc")
xlabel('Frequency [Hz]')
ylabel('|H| [m/s*N]')
grid on
axis tight
legend('sc-sc','oc-sc');
