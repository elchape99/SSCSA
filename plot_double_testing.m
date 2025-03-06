close all
pippo = 2400;


FRF_sc_sc_double_piezo_fin = load ("FRF/FRF_sc_sc_double_piezo_final.mat");
FRF_rl_rl_temp = load ("FRF/FRF_rl_rl_double_piezo_interm.mat");
figure
semilogy (freq(:,1:pippo), FRF_sc_sc_double_piezo_fin.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl(1:pippo,:))
hold on
semilogy (freq(:,1:pippo), FRF_rl_rl_temp.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl(1:pippo,:))

% H_rl_rl_test_1_1 = 1i .* w .* double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], ...
%     [beam.Cp.C21, beam.Cp.C22], L1_opt_new, L2_opt_new, 0.5 * R1_opt_new, R2_opt_new, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);
% 
% H_rl_rl_test_1_2 = 1i .* w .* double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], ...
%     [beam.Cp.C21, beam.Cp.C22], L1_opt_new, L2_opt_new, 1.5 * R1_opt_new, R2_opt_new, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);
% 
% 
% H_rl_rl_test_2_1 = 1i .* w .* double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], ...
%     [beam.Cp.C21, beam.Cp.C22], 0.9 * L1_opt_new, L2_opt_new, R1_opt_new, R2_opt_new, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);
% 
% H_rl_rl_test_2_2 = 1i .* w .* double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], ...
%     [beam.Cp.C21, beam.Cp.C22], 1.1 * L1_opt_new, L2_opt_new, R1_opt_new, R2_opt_new, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);
% 
% H_rl_rl_test_3 = 1i .* w .* double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], ...
%     [beam.Cp.C21, beam.Cp.C22], 0.8 * L1_opt_new, L2_opt_new, 0.58 * R1_opt_new, R2_opt_new, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);
% 
% H_rl_rl_test_4 = 1i .* w .* double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], ...
%     [beam.Cp.C21, beam.Cp.C22], L1_opt_new, L2_opt_new, 0.58 * R1_opt_new, R2_opt_new, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);
% 
% H_rl_rl_test_5 = 1i .* w .* double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], ...
%     [beam.Cp.C21, beam.Cp.C22], L1_opt_new, L2_opt_new, R1_opt_new, R2_opt_new, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);
% 
% hold on 
% % semilogy (freq(:,1:pippo), abs(sum(H_rl_rl_test_1_1(1:pippo,:), 2)))
% % semilogy (freq(:,1:pippo), abs(sum(H_rl_rl_test_1_2(1:pippo,:), 2)))
% % 
% % semilogy (freq(:,1:pippo), abs(sum(H_rl_rl_test_2_1(1:pippo,:), 2)))
% % semilogy (freq(:,1:pippo), abs(sum(H_rl_rl_test_2_2(1:pippo,:), 2)))
% 
% semilogy (freq(:,1:pippo), abs(sum(H_rl_rl_test_3(1:pippo,:), 2)))
% semilogy (freq(:,1:pippo), abs(sum(H_rl_rl_test_4(1:pippo,:), 2)))
% legend("sc sc", "rl rl", "mix", "L corr")
% 
% figure 
% plot (freq(1:pippo), 20 * log10(FRF_sc_sc_double_piezo_fin.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl(1:pippo)))
% hold on 
% plot (freq(1:pippo), 20 * log10(abs(sum(H_sc_sc_fitted(1:pippo, :), 2))))
% 
% plot (freq(1:pippo), 20 * log10(abs(sum(H_rl_rl_test_4(1:pippo, :),2))))
% 
% plot (freq(1:pippo), 20 * log10(abs(sum(H_rl_rl_test_5(1:pippo, :),2))))


H_rl_rl_test_6 = 1i .* w .* double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], ...
    [beam.Cp.C21, beam.Cp.C22], L1_opt_new, 1.1*L2_opt_new, R1_opt_new, R2_opt_new, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);

H_rl_rl_test_7 = 1i .* w .* double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], ...
    [beam.Cp.C21, beam.Cp.C22], L1_opt_new, 0.95*L2_opt_new, R1_opt_new, R2_opt_new, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);

H_rl_rl_test_8 = 1i .* w .* double_piezo_reson_FRF (w, w_i, w_cap, csi_i, [beam.Cp.C11, beam.Cp.C12], ...
    [beam.Cp.C21, beam.Cp.C22], L1_opt_new, L2_opt_new, R1_opt_new, R2_opt_new, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);

semilogy (freq(:,1:pippo), abs(sum(H_rl_rl_test_6(1:pippo,:), 2)))
semilogy (freq(:,1:pippo), abs(sum(H_rl_rl_test_7(1:pippo,:), 2)))
legend("sc sc", "rl rl", "1.1 L", "0.9 L")