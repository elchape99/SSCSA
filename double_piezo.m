clc
clear all
%% Script developed for use two piezo

csi_i = 0; % natural damping of the system
% with different damping change the peak, but the intersection point with
% the short circuit response does not change
csi1 = 1;  
csi2 = 1;
w = 1:6000;
% w_i: freq naturali
% w_1: freq ottime piezo 1
% w_2: freq ottime piezo 2

% response function with two piezo with resonant
H_rl_rl = (w^6 + ...
    2* 1i * w.^5 * (csi1 * w_2 + csi2 * w_2 + csi_i * w_i)...
    + w.^4 * (w_2.^2 + w_1.^2 + 4 * csi1 * csi2 * w_1 * w_2 + 4 * csi_i * w_i (csi_1 * w_1 + csi2 * w_2) + w_i.^2) ...
    - 2 * 1i * w.^3 * (csi1 * w_1 * w_2.^2 + csi2 * w_2 * w_1.^2 + csi_i * w_i * (w2.^2 + w1.^2 + 4 * csi1 * csi2 * w_1 * w_2 + w_i.^2 * (csi1* w_1 + csi2 * w_2)))...
    - w.^2 * (w_1.^2 * w_2.^2 + 4 * csi_i * w_i * (csi1 * w_1 * w_2.^2 + csi2 * w_2 * w_1.^2) + w_i * (w_1.^2 + w_2.^2 + 4 * csi1*csi2 * w_1 * w_2 + k1 * w_1.^2 + k2 * w_2.^2))...
    + 2 * 1i * w_i * w * (csi_i * w_i * w_1.^2 * w_2.^2 + csi1 * w_1 * w_2.^2 + csi2 * w_2 * w_1.^2 - csi2 * k1.^2 * w_1.^2 * w_2 - csi2 * k2.^2 * w_1 * w_2.^2) ...
    + w_i.^2 * (w_1.^2 * w_2.^2 + k1.^2 * w_1.^2 * w_2.^2 + k1.^2 * w_1.^2 * w_2.^2 + k2.^2 * w_1^2 * w_2.^2))...
    /((-w.^2 + 2 * 1i * csi1 * w_1 * w + w_1.^2) * (-w.^1 * 2 * 1i * csi2 * w_2 * w + w_2.^2));
      
H_sc_sc = 1/(w_i^2)