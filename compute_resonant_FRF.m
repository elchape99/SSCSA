function [reson_FRF] = compute_resonant_FRF(w, w_i, w_i_oc, csi_k, w_e, csi_e)
        
        reson_FRF_num = -w.^2 + w_e.^2 + 2j *csi_e .* w_e .* w;
        reson_FRF_den_1 = w.^4 + w_i.^2 .* w_e.^2;
        reson_FRF_den_2 = - w.^2 .* (w_e.^2 + 4 * csi_k .* csi_e .* w_i .* w_e + w_i_oc.^2);
        reson_FRF_den_3_a = 2 * csi_e .* w_e .* (w_i_oc.^2 - w.^2);
        reson_FRF_den_3_b = 2 * csi_k .* w_i .* (w_e.^2 - w.^2);
        reson_FRF_den_3 = 1j * w .* (reson_FRF_den_3_a + reson_FRF_den_3_b);
        reson_FRF_den = reson_FRF_den_1 + reson_FRF_den_2 + reson_FRF_den_3;

        reson_FRF = reson_FRF_num ./ reson_FRF_den;
end

