clear 
close all
clc
addpath("Functions\");

%% RESISTIVE CIRCUIT
% load beam data
beam_defm
nmodes=3;

% natural pulse of the system, no action of the piezo 
w_i = 2*pi*beam.nf.sc_sc; 

% ki = sqr ((w_i_oc ^ 2 - w_i.^2)/w_i.^2)
beam.k.k1 = sqrt ((beam.nf.oc_sc .^2 - beam.nf.sc_sc.^2) ./ beam.nf.sc_sc .^2);
beam.k.k2 = sqrt ((beam.nf.sc_oc .^2 - beam.nf.sc_sc.^2) ./ beam.nf.sc_sc .^2);

% find tau optimal, piezo 1, mode 1, resistant otpimization
beam.wf.wf1 = sqrt (w_i.^2 .* (1 + beam.k.k1.^2 ./ 2));
beam.tau.tau1 = 1 ./ beam.wf.wf1;
beam.R.R1 = beam.tau.tau1 ./ beam.Cp.C11; 

% find tau optimal, piezo 2, mode 2, resonant optimization
beam.wf.wf2 = sqrt (w_i.^2 .* (1 + beam.k.k2.^2 ./ 2));
beam.tau.tau2 = 1 ./ beam.wf.wf2;
beam.R.R2 = beam.tau.tau2 ./ beam.Cp.C12; 


FRF.sc_sc = load ("FRF/FRF_sc_sc.mat");
FRF.oc_sc = load ("FRF/FRF_oc_sc.mat");
FRF.sc_oc = load ("FRF/FRF_sc_oc.mat");
FRF.R_sc = load("FRF/FRF_R_sc.mat");
FRF.sc_sc_new = load ("FRF/FRF_sc_sc_new");
FRF.sc_RL = load ("FRF/FRF_sc_RL");

max_sc_sc = max(FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl);
max_R_sc = max(FRF.R_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl);


precentage = ((max_sc_sc - max_R_sc)/ max_sc_sc) * 100;
percentage2 = ((0.194237 - 0.0167)/0.194237) * 100;

% x line for the plot
x = linspace(1,250, length(FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl(1:3000)));

% % plot the total FRF, all processing for mode 1 (sc, oc, only pezo 1)
% figure
% semilogy (x, FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl, "LineWidth",0.6);
% hold on
% semilogy(x, FRF.oc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl, "LineWidth",0.6);
% hold on
% title ("FRF-sc-sc - FRF-oc-sc")
% xlabel('Frequency [Hz]')
% ylabel('|H| [m/s*N]')
% grid on
% axis tight
% legend('sc-sc','oc-sc');
% 
% % First piezo optimized with Resistance
% figure
% semilogy (x, FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl, "LineWidth",0.3);
% hold on
% semilogy(x, FRF.oc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl, "LineWidth",0.3);
% hold on
% semilogy(x, FRF.R_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl, "LineWidth",1, "Color","g");
% title ("FRF bottom piezo with Resistance Optimization")
% xlabel('Frequency [Hz]')
% ylabel('|H| [m/s*N]')
% grid on
% axis tight
% legend('sc-sc','oc-sc', 'R-sc');
size(FRF.sc_sc_new.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl)
% Second piezo optimized with R + L 
figure
semilogy (x, FRF.sc_sc_new.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl(1:3000), "LineWidth",0.6);
hold on
semilogy(x, FRF.sc_RL.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl(1:3000), "LineWidth",0.6);
hold on
title ("FRF upper piezo with R + L  Optimization")
xlabel('Frequency [Hz]')
ylabel('|H| [m/s*N]')
grid on
axis tight
legend('sc-sc','sc-RL');





