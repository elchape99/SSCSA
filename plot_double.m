%% Script for plot all the graph needed for double piezo 
close all
% values from double_piezo.m script
% Run this scropt after runned double_piezo.m script

%% === Plot the fitted response ===
figure(1);
subplot(2, 1, 1);
plot(freq, FRF_sc_sc, 'r--'); 
hold on;
plot(freq, H_sc_sc_fitted, 'b-');
xlabel('Frequency (Hz)'); 
ylabel('|H(ω)|');
axis tight
grid on;
title('Curves fitting on decimal scale');
legend('FRF-sc-sc', 'H-sc-sc-fitted');


subplot(2, 1, 2);
semilogy (freq, (FRF_sc_sc));
hold on
semilogy(freq, (H_sc_sc_fitted));
xlabel('Frequency [Hz]')
ylabel('|H| [m/(s*N)]')
axis tight
grid on
title ('Curves fitting on logaritmic scale')
legend('FRF-sc-sc','H-sc-sc-fitted');

%% Plot the three finded FRF/H_rl_rl without optimization
figure(2);
semilogy (freq, (FRF_rl_rl), "LineWidth",0.6);
hold on
semilogy(freq, abs(sum(H_rl_rl,2)));
semilogy(freq, abs(H_rl_rl_fitted));
xlabel('Frequency [Hz]')
ylabel('|H| [m/(s*N)]')
axis tight
grid on
title('Comparison between all FRF/H rl-rl without optimization');
legend ('FRF-rl-rl', 'H-rl-rl', 'H-rl-rl-fitted');

%% plotting the H (normal vs fitted)
figure(3);
subplot(3,1,1)
semilogy(freq, abs(sum(H_sc_sc, 2)))
hold on
semilogy(freq, abs(sum(H_rl_rl,2)));
%semilogy(w, abs(H_rl_rl), w, abs(sum(H_sc_sc,2)));
%legend("H-rl-rl", "H-rl-rl-2", "H-sc-sc")
xlabel('Frequency [Hz]')
ylabel('|H| [m/(s*N)]')
grid on;
axis tight
title('analitic plot of H-sc-sc and H-rl-rl')
legend("H-sc-sc", "H-rl-rl")

% plot the total FRF, all processing for mode 1
subplot(3,1,2)
semilogy (freq, H_sc_sc_fitted, "LineWidth",0.6);
hold on
semilogy(freq, H_rl_rl_fitted, "LineWidth",0.6);
xlabel('Frequency [Hz]')
ylabel('|H| [m/(s*N)]')
axis tight
grid on
title ("all measured value FRF-sc-sc - FRF-rl-rl")
legend('FRF-sc-sc','FRF-rl-rl');

subplot(3, 1, 3)
semilogy (freq, H_sc_sc_fitted, "LineWidth",0.6);
hold on
semilogy(freq, H_rl_rl_fitted, "LineWidth",0.6);
semilogy(freq, abs(sum(H_sc_sc, 2)))
semilogy(freq, abs(sum(H_rl_rl, 2)))
xlabel('Frequency [Hz]')
ylabel('|H| [m/(s*N)]')
axis tight
grid on
title ("all measured value FRF-sc-sc - FRF-rl-rl")
legend('H-sc-sc-fitted','H-rl-rl-fitted','H-sc-sc','H-rl-rl');

%% Optimization L plot
% Visualizzazione della funzione con i parametri ottimali
figure(4);
subplot (2, 1, 1);
plot(w, best_FRF_values, 'b');
hold on;
plot(best_locs, best_pks, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
xlabel('Pulse (w = 2*pi*f)');
ylabel('FRF Magnitude');
axis tight
grid on;
title('FRF module with optimal parameters');
legend('FRF', 'First 4 peaks');

subplot(2, 1, 2)
semilogy(freq, (abs(1i .* w .* best_FRF_values)));
hold on
xline (best_locs / (2 * pi))
semilogy(freq, abs(H_sc_sc_fitted))
semilogy(freq, abs(H_rl_rl_damped))
semilogy (freq, abs(H_rl_rl_damped_opt))
xlabel('Frequency (Hz)');
ylabel('|H| [m/(s*N)]');
axis tight
grid on;
title('FRF optimal in logaritmic scale with optimal values');
legend('H-rl-rl L opt', 'H-sc-sc-fitted', 'H-rl-rl-damped', 'H-rl-rl-opt-damped');

figure
semilogy(freq, abs(H_sc_sc_fitted))
hold on
semilogy (freq, abs(H_rl_rl_damped_opt))
grid on
legend ('H-sc-sc', 'H-rl-rl-optimal')

figure
H_sc_sc_unnatdamp = phi_opt .* 1i .* w ./ (w_i.^2 - w.^2);
H_rl_rl_unnatdamp = double_piezo_reson_FRF(w, w_i, w_cap, csi_i_optimization, ...
                              [beam.Cp.C11, beam.Cp.C12], [beam.Cp.C21, beam.Cp.C22], ...
                              L1_opt_new, L2_opt_new, R1_opt_new, R2_opt_new, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);
H_rl_rl_unnatdamp = phi_opt .* 1i .* w .* H_rl_rl_unnatdamp;

semilogy(freq, abs((H_rl_rl_fitted)))
hold on
semilogy (freq, abs(sum(H_rl_rl_damped_opt, 2)))
grid on
legend ('H-rl-rl', 'H-rl-rl-optimal')