clc
close all
clear 
%% Script developed for use two piezo
% name convenction:
% - H: analitic function, computed with formulas
% - FRF: experimental measured valuesù
% - w_cap = w_oc_oc
% - w_i: freq naturali
% - w_1: freq ottime piezo 1
% - w_2: freq ottime piezo 2


% data of the bar 
beam_defm;

% load experimantal data
FRF.sc_sc = load ("FRF/FRF_sc_sc.mat");
FRF_rl_rl_double_piezo_first = load ("FRF/FRF_rl_rl_double_piezo_first");

% rename correctly the variable
FRF_sc_sc = FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl;
FRF_rl_rl = FRF_rl_rl_double_piezo_first.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl;

csi_i = beam.xi.sc_sc(1:2); % natural damping of the system
% with different damping change the peak, but the intersection point with
% the short circuit response does not change
csi1 = 0.005;  
csi2 = 0.005;  
k = [0.24, 0.18];
k1= k(1);
k2= k(2);

% freq used for the plot
freq = linspace(1,500, length(FRF_sc_sc)); 

% trasform freq in pulse
w = freq' * 2 * pi;

% w_i = beam.nf.sc_sc(1:2); % w_sc_sc
% w_i = [21.1297, 116.538];
% w_cap = beam.nf.oc_oc(1:2); % w_oc_oc
w_i = 2*pi* beam.nf.sc_sc(1:2); % w_sc_sc

w_i = [21.1297*2*pi, 116.538 * 2*pi]; % trovate in matlab
w_cap = 2*pi* beam.nf.oc_oc(1:2); % w_oc_oc

%optimal values if consider two modes uncorrelated 
w_e_1_opt = (w_cap(1));
w_e_2_opt = (w_cap(2));

csi_e_1_opt = sqrt(3)/2 * sqrt ((w_cap(1)^2 - w_i(1)^2) / (w_cap(1)^2 + w_i(1)^2));
csi_e_2_opt = sqrt(3)/2 * sqrt ((w_cap(2)^2 - w_i(2)^2) / (w_cap(2)^2 + w_i(2)^2));

% computed with w_e_1 is a pulse not a period
L1_opt = 1 ./ (w_e_1_opt^2 * (beam.Cp.C11));
L2_opt = 1 ./ (w_e_2_opt^2 * (beam.Cp.C22));

% L in freq
L1_opt = L1_opt/(2*pi)^2;
L2_opt = L2_opt/(2*pi)^2;


% R1_opt_1 = 2 * csi_e_1_opt * sqrt (L1_opt / beam.Cp.C11);
% R2_opt_1 = 2 * csi_e_2_opt * sqrt (L2_opt / beam.Cp.C22);
R1_opt = (2 * csi_e_1_opt) / (beam.Cp.C11 * w_e_1_opt);
R2_opt = (2 * csi_e_2_opt) / (beam.Cp.C22 * w_e_2_opt);

% R in freq
R1_opt = R1_opt * 2*pi;
R2_opt = R2_opt * 2*pi;

% analitic H_sc_sc
H_sc_sc = 1./(w_i.^2 + 1i.*2.*csi_i .*w_i.*w - w.^2);
% experimental measure of velocity force : moltiplication of all H(w)
% analitic by jw
H_sc_sc = 1i * w .* H_sc_sc;

% analitic H_rl_rl
H_rl_rl = double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], [beam.Cp.C21, beam.Cp.C22], L1_opt, L2_opt, R1_opt, R2_opt, beam.k.k1(1:2), beam.k.k2(1:2));      
% experimental measure of velocity7force : moltiplication of all H(w)
% analitic by jw
H_rl_rl = 1i * w .* H_rl_rl;

%% Curve fitting

% === Funzione modello della FRF ===
H_model = @(phi_squared, x) abs(phi_squared(1) .* H_sc_sc(:,1) + ...
                                phi_squared(2) .* H_sc_sc(:,2));

% === Curve fitting sulla magnitudo della FRF ===
phi_init = [0.1, 0.1]; % Stima iniziale dei pesi modali
phi_opt = lsqcurvefit(H_model, phi_init, double(freq), double(abs(FRF_sc_sc)));

H_sc_sc_fitted = abs(phi_opt(1) * H_sc_sc(:,1) + ...
                     phi_opt(2) * H_sc_sc(:,2));

% === Risultati ===
disp('Mode shapes relativi nel punto di misura:');
disp(sqrt(phi_opt));

% === Plot del risultato ===
figure;
subplot(2, 1, 1);
plot(freq, FRF_sc_sc, 'r--', 'DisplayName', 'FRF sperimentale'); 
hold on;
plot(freq, H_model(phi_opt, freq), 'b-', 'LineWidth', 2, 'DisplayName', 'FRF fittata');
xlabel('Frequenza (Hz)'); ylabel('|H(ω)|');
legend('FRF-sc-sc', 'H_model');
grid on;
title('Curve fitting della FRF');

subplot(2, 1, 2);
semilogy (freq, (FRF_sc_sc), "LineWidth",0.6);
hold on
semilogy(freq, (H_sc_sc_fitted));
%semilogy(w, abs(sum(H_sc_sc,2)));
title ("H-sc-sc-fitted - FRF-sc-sc", 'H-sc-sc')
xlabel('Frequency [Hz]')
ylabel('|H| [m/s*N]')
grid on
axis tight
legend('FRF-sc-sc','H-sc-sc-fitted');

H_rl_rl_fitted = abs(phi_opt(1) * H_rl_rl(:,1) + ...
                     phi_opt(2) * H_rl_rl(:,2));

figure;
semilogy (freq, (FRF_rl_rl), "LineWidth",0.6);
hold on
semilogy(freq, (H_rl_rl_fitted));
semilogy(freq, abs(sum(H_rl_rl,2)));
title ("FRF-rl-rl", 'H-rl-rl-fitted')
xlabel('Frequency [Hz]')
ylabel('|H| [m/s*N]')
grid on
axis tight
legend('FRF-rl-rl','H-rl-rl-fitted');

% plotting the H (analitic frequecies respone function)
figure;
subplot(2,1,1)
semilogy(w, abs(sum(H_sc_sc, 2)), w, abs(sum(H_rl_rl,2)));
legend("H-sc-sc", "H-rl-rl")
title('analitic plot of H-sc-sc and H-rl-rl')
%semilogy(w, abs(H_rl_rl), w, abs(sum(H_sc_sc,2)));
%legend("H-rl-rl", "H-rl-rl-2", "H-sc-sc")
xlabel('Frequency [Hz]')
ylabel('|H| [m/s*N]')
grid on;
axis tight

% plot the total FRF, all processing for mode 1 (sc, oc, only pezo 1)
subplot(2,1,2)
semilogy (freq, FRF_sc_sc, "LineWidth",0.6);
hold on
semilogy(freq, FRF_rl_rl, "LineWidth",0.6);
hold on
title ("all measured value FRF-sc-sc - FRF-rl-rl")
legend('FRF-sc-sc','FRF-rl-rl');
xlabel('Frequency [Hz]')
ylabel('|H| [m/s*N]')
grid on
axis tight


% % overlapping of the two FRFsc-sc (analitic vs measured)
% figure;
% subplot(2,1,1)
% semilogy (x, FRF_sc_sc, "LineWidth",0.6);
% hold on
% semilogy(w, abs(sum(H_sc_sc,2)));
% hold on
% title ("FRF-sc-sc - H-sc-sc")
% xlabel('Frequency [Hz]')
% ylabel('|H| [m/s*N]')
% grid on
% axis tight
% legend('FRF-sc-sc','H-sc-sc');

% find the peack and normalization
[magnitudo_p_m, position_p_m] = findpeaks(FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl,'MinPeakHeight',0.1);
[magnitudo_p_a, position_p_a] = findpeaks(abs(sum(H_sc_sc(:,1),2)));
peak = max(H_sc_sc(:,1));

H_sc_sc_norm = H_sc_sc ./ magnitudo_p_a(1);
FRF.sc_sc_norm = FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl ./ magnitudo_p_m(1);





%%
% 
% 
% 
% 
% 
% % optimization of the analitic peak
% 
% % % [a; b] * [h_sc_sc_1 (w1), h_sc_sc_2(w1); h_sc_sc_2(w2), h_sc_sc_2(w2)] =
% % % [p_m_1; p_m_2]
% % peak_idx = abs(w-w_i) < 0.05;
% % H_cand_11 = H_sc_sc ( peak_idx(:, 1), 1);
% % H_cand_12 = H_sc_sc(peak_idx (:, 1), 2);
% % H_cand_21 = H_sc_sc(peak_idx (:, 2), 1);
% % H_cand_22 = H_sc_sc(peak_idx (:, 2), 2);
% % 
% % PHI_squared = ([H_cand_11(1), H_cand_12(1); H_cand_21(1), H_cand_22(1)])\magnitudo_p_m;
% % PHI_squared = magnitudo_p_m ./ magnitudo_p_a;
% 
% 
% peak_idx = abs(w-w_i) < 0.05 * 2* pi;
% peak_idx (find(peak_idx(:, 1), 1) + 1 : end, 1) = false;
% peak_idx (find(peak_idx(:, 2), 1) + 1 : end, 2) = false;
% pippo1 = H_sc_sc_norm(peak_idx);
% pippo2 = abs(FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl(peak_idx(:, 1)' | peak_idx(:, 2)'));
% 
% % sequenz: abs-sum-abs
% H_diff = @(PHI_sq) abs(sum(abs(FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl(peak_idx(:, 1)' | peak_idx(:, 2)')) - abs(PHI_sq .* H_sc_sc(peak_idx))));
% H_diff_norm = @(PHI_sq) abs(sum(abs(FRF.sc_sc_norm(peak_idx(:, 1)' | peak_idx(:, 2)')) - abs(PHI_sq .* H_sc_sc_norm (peak_idx))));
% 
% PHI_squared = fminsearch (H_diff, [1; 1]);
% PHI_squared_norm = fminsearch (H_diff_norm, [1; 1]);
% 
% % compute the new H with the weight
% H_sc_sc_weighted = PHI_squared(1) .* H_sc_sc(:,1) + PHI_squared(2) .* H_sc_sc(:,2); 
% H_sc_sc_weighted_norm = PHI_squared_norm(1) .* H_sc_sc_norm(:,1) + PHI_squared_norm(2) .* H_sc_sc_norm(:,2); 
% H_rl_rl_weighted = PHI_squared(1) * H_rl_rl(:,1) + PHI_squared(2) * H_rl_rl(:,2);
% 
% 
% figure(4)
% semilogy (x, abs(sum(H_sc_sc_weighted, 2)), "LineWidth",0.6);
% hold on
% semilogy(w, abs(FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl));
% title ("H-sc-sc-weighted - FRF-sc-sc")
% xlabel('Frequency [Hz]')
% ylabel('|H| [m/s*N]')
% grid on
% axis tight
% legend('H-sc-sc-weighted', 'FRF-sc-sc');
% 
% figure(5)
% semilogy (x, abs(sum(H_sc_sc_weighted_norm, 2)), "LineWidth",0.6);
% hold on
% semilogy(w, abs(FRF.sc_sc_norm));
% title ("H-sc-sc-weighted-norm - FRF-sc-sc_norm")
% xlabel('Frequency [Hz]')
% ylabel('|H| [m/s*N]')
% grid on
% axis tight
% legend('H-sc-sc-weighted-norm', 'FRF-sc-sc-norm');
% 
% figure(6)
% semilogy (x, abs(H_rl_rl_weighted), "LineWidth",0.6);
% hold on
% semilogy(w, abs(sum(H_rl_rl, 2)));
% title ("H-rl-rl-weighted - H-rl-rl")
% xlabel('Frequency [Hz]')
% ylabel('|H| [m/s*N]')
% grid on
% axis tight
% legend('H-rl-rl-weighted', 'H-rl-rl');
% 
% figure(7)
% semilogy (x, abs(sum(H_rl_rl,2)), "LineWidth",0.6);
% hold on
% semilogy(w, abs(FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl));
% title ("H-rl-rl - FRF-rl-rl")
% xlabel('Frequency [Hz]')
% ylabel('|H| [m/s*N]')
% grid on
% axis tight
% legend('H-rl-rl-weighted', 'FRF-rl-rl');
% 
% [mode_shapes, FRF_weighted] = extract_modeshape_weights(FRF.sc_sc.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl, x, H_sc_sc(:,1), x);
% 
% figure(8)
% semilogy (x, abs(FRF_weighted), "LineWidth",0.6);
% hold on
% semilogy(w, abs(H_sc_sc));
% title ("H-rl-rl - FRF-rl-rl")
% xlabel('Frequency [Hz]')
% ylabel('|H| [m/s*N]')
% grid on
% axis tight
% legend('FRF-w-1','FRF-w-2','H-sc-sc-1', 'H-sc-sc-2');
% 
% 
