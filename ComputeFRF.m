function [Hi, H] = ComputeFRF(phi, nf, omega, csi, m, f)
% Function for comupute Frequency response function
% Input:
% - phi: mode shape matrix 
% - nf: natural frequences of the system
% - omega: range of system excitation frequencies
% - csi: adimensional damping ratio 
% - m: measurement point
% - f = forcing point
%
% Output
% - H: frequency response function asscotated at N mode
n_mode = size(phi, 2);

% compute the natural frequences (hz -> rad/s)
omega_n = 2*pi*nf;
Hi = zeros(length(omega), n_mode); % the lenght of H is equal to the lenght of 
H  = zeros(length(omega), 1);

for r = 1:n_mode % for cycle conidering r-th mode
    Hi(:,r) = phi(m,r)*phi(f,r)*((1./(-omega.^2 + 2*1i*omega_n(r).*csi(r).*omega + omega_n(r).^2))); 
    H = H + Hi(:, r);
end