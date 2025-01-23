function [cont_FRF] = compute_resistive_FRF (w, w_i, w_i_oc, csi_k, tau)
    %w has number of modes rows
    % cont_FRF = zeros(size(w));
    % H_cont(jw) = 1 + jwtau/ (wi^2 - w^2*(1+2taucsiiwi) + jw(tauwoci^2 +
        % 2wicsii - tau w^2)
    %denomintor first two terms
    cont_FRF_den = j * w .* (tau .* w_i_oc.^2 + 2 * w_i .* csi_k - tau .* w.^2);
    cont_FRF_den = cont_FRF_den + w_i.^2 - w.^2 .* (1 + 2 * csi_k .* tau .* w_i);
    % cont_FRF_den = cont_FRF_den + 
    cont_FRF = (1 + j * w .* tau) ./ cont_FRF_den;
end

