%%
clear
clc
close all
load Data.mat
%%
% csi = [csi1; csi2; csi3];
% nfi_oc = [nf1_oc; nf2_oc; nf3_oc];
% nfi_sc = [nf1_sc; nf2_sc; nf3_sc];

%%


%%
%sensing at x = L/6, foce at x = 5/6 L
s_idx = round (length (PHI) / 6);
f_idx = round (5/6 * length (PHI));

% w = 1:0.0001:1000;
% w = 1:1000;
w = 1:0.01:1000;

mode_num = 4;
modes = 1:4;
single_FRF_i = zeros (mode_num, length (w));
single_FRF_i_oc = zeros (mode_num, length (w));
total_FRF = zeros (1, length (w));
total_FRF_oc = zeros(1, length(w));
% ki = zeros (1, mode_num);
col_w = ones (mode_num, 1) .* w;
w_i = zeros (1, mode_num);
w_i_oc = zeros(1, mode_num);
csi_i = zeros (1, mode_num);

csi_i = csi1;
w_i = nf1_sc;
w_i_oc = nf1_oc;

ki = sqrt ((w_i_oc.^2 - w_i.^2) / w_i.^2);

single_FRF_i (:, :)  = 1./(w_i'.^2 - col_w.^2 + 1i * 2 * csi_i' .* w_i' .* col_w);
single_FRF_i_oc (:, :) = 1./(w_i_oc'.^2 - col_w.^2 + 1i * 2 * csi_i' .* w_i_oc' .* col_w);


total_FRF = sum (PHI(f_idx, modes)' .* PHI(s_idx, modes)' .* single_FRF_i (:, :));
total_FRF_oc (:) = sum (PHI(f_idx, modes)' .* PHI(s_idx, modes)' .* single_FRF_i_oc (:, :));

%compute optimal resistance
%w_f^2 = w_i^2 * (1 + ki^2 / 2)
%tau = R*Cpi
% tau_opt = 1/w_f

w_f = sqrt (w_i.^2 .* (1 + ki.^2 ./ 2));
tau_opt = 1 ./ w_f;
R_opt = tau_opt ./ cpi_tot;

w_e_opt = w_i_oc;
csi_e_opt = sqrt(3)/2 * sqrt ((w_i_oc.^2 - w_i.^2) ./ (w_i_oc.^2 + w_i.^2))';

L_opt (:) = 1 ./ (w_e_opt.^2 .* cpi_tot);
%optimal resistance for series
R_reson_opt = 2 * csi_e_opt ./ (w_e_opt .* cpi_tot);
%optimal resistance for parallel
% R_reson_opt = 1 ./ (2 * csi_e_opt .* w_e_opt .* cpi_tot);

cont_mode = 1;
resist_FRF_set = compute_resistive_FRF (col_w, w_i', w_i_oc', csi_i', tau_opt (cont_mode));
resist_FRF = sum (PHI(f_idx, modes)' .* PHI(s_idx, modes)' .* resist_FRF_set);

reson_FRF_set = compute_resonant_FRF (col_w, w_i', w_i_oc', csi_i', w_e_opt (cont_mode), csi_e_opt (cont_mode));
        
reson_FRF = sum (PHI(f_idx, modes)' .* PHI(s_idx, modes)' .* reson_FRF_set);


%point 1 plotting
figure
subplot (2, 1, 1)
semilogy(w, abs(single_FRF_i))
hold on
semilogy(w, abs(total_FRF));
legend ('mode 1', 'mode 2', 'mode 3', 'mode 4', 'total')
subplot (2, 1, 2)
plot (w, angle (single_FRF_i))
hold on
plot (w, angle (total_FRF))
legend ('mode 1', 'mode 2', 'mode 3', 'mode 4', 'total')

% point 2-5
%point 2 plotting
figure
subplot (2, 1, 1)
semilogy (w, abs(total_FRF), w, abs(total_FRF_oc));
legend ('response short circuit', 'response open circuit')
subplot (2, 1, 2)
plot (w, angle (total_FRF), w, angle (total_FRF_oc));
legend ('response short circuit', 'response open circuit')

%point 3
% disp('electro-mechanical coupling coefficients: ')
% disp (ki(k, :))
% w_i = nfi_sc (k, mode);
% w_i_oc = nfi_oc (k, mode);
% ki(k, mode) = sqrt ((w_i_oc^2 - w_i^2)/w_i^2);


%point 4  plotting
figure ('Name', "FRF oc and sc compared with resistive and resonant shunt")
% loglog(w, abs(total_FRF), w, abs(total_FRF_oc));
subplot (2, 1, 1)
semilogy(w, abs(total_FRF), w, abs(total_FRF_oc));
hold on
semilogy (w, abs (resist_FRF))
semilogy (w, abs (reson_FRF))
% axis ([0.85 * w_i - 20, 1.15 * w_i + 20, 0, 1.2])
% axis ([w_i - 20, w_i + 20, 0, 1.2])
% xline (w_f(k, M))
hold off
legend ('sc', 'oc', 'opt resist', 'opt reson')
subplot (2, 1, 2)
plot (w, angle(total_FRF), w, angle(total_FRF_oc));
hold on
plot (w, angle (resist_FRF))
plot (w, angle (reson_FRF))
hold off
legend ('sc', 'oc', 'opt resist', 'opt reson')
