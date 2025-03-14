function [H_rl_rl] = double_piezo_reson_FRF(w, w_i, w_cap, csi_i, C1i, C2i, L1, L2, R1, R2, k1, k2, PHI_sqr)

w_1 = 1./sqrt(C1i .* L1);
w_2 = 1./sqrt(C2i .* L2);
csi1 = R1./ 2 .* sqrt (C1i ./ L1);
csi2 = R2./ 2 .* sqrt (C2i ./ L2);

H_rl_rl = (PHI_sqr .* (-w.^2 + 2 .* 1i .* csi1 .* w_1 .* w + w_1.^2) .* (-w.^2 + 2 .* 1i .* csi2 .* w_2 .* w + w_2.^2)) ...
    ./ ...
    (-w.^6 + ...
    2.* 1i .* w.^5 .* (csi1 .* w_1 + csi2 .* w_2 + csi_i .* w_i) ...
    + w.^4 .* (w_2.^2 + w_1.^2 + 4 .* csi1 .* csi2 .* w_1 .* w_2 + 4 .* csi_i .* w_i .* (csi1 .* w_1 + csi2 .* w_2) + w_cap.^2) ...
    - 2 .* 1i .* w.^3 .* (csi1 .* w_1 .* w_2.^2 + csi2 .* w_2 .* w_1.^2 + csi_i .* w_i .* (w_2.^2 + w_1.^2 + 4 .* csi1 .* csi2 .* w_1 .* w_2) + w_cap.^2 .* (csi1.* w_1 + csi2 .* w_2))...
    - w.^2 .* (w_1.^2 .* w_2.^2 + 4 .* csi_i .* w_i .* (csi1 .* w_1 .* w_2.^2 + csi2 .* w_2 .* w_1.^2) + w_cap.^2 .* (w_1.^2 + w_2.^2 + 4 .* csi1.*csi2 .* w_1 .* w_2) - w_i.^2 .* (k1.^2 .* w_1.^2 + k2.^2 .* w_2.^2)) ...
    + 2 .* 1i .* w .* (csi_i .* w_i .* w_1.^2 .* w_2.^2 + w_cap.^2 .* (csi1 .* w_1 .* w_2.^2 + csi2 .* w_2 .* w_1.^2) - w_i.^2 .* (csi2 .* k1.^2 .* w_1.^2 .* w_2 + csi1 .* k2.^2 .* w_1 .* w_2.^2)) ...
    + w_cap.^2 .* (w_1.^2 .* w_2.^2) - w_i.^2 .* (k1.^2 .* w_1.^2 .* w_2.^2 + k2.^2 .* w_1.^2 .* w_2.^2));
end