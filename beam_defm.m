% bar definition
% - first piezo is the one on the bottom
% - second piezo id the one on the top
clc
close all

% config with sc_sc
% da software
beam.nf.sc_sc(1) = 20.122; %[Hz]
beam.nf.sc_sc(2) = 115.74; %[Hz]
beam.nf.sc_sc(3) = 314.24; %[Hz]

% sc_sc_double_piezo
beam.nf.sc_sc_double_piezo(1) = 20.092; %[Hz]
beam.nf.sc_sc_double_piezo(2) = 115.6; %[Hz]


% da software
beam.xi.sc_sc(1) = 0.010033;
beam.xi.sc_sc(2) = 0.0047168;
beam.xi.sc_sc(3) = 0.0028293;

% config with oc_sc (piezo 1)
beam.nf.oc_sc(1) = 20.723; %[Hz]
beam.nf.oc_sc(2) = 115.94; %[Hz]
beam.nf.oc_sc(3) = 314.81; %[Hz]

beam.xi.oc_sc(1) = 0.0086232;
beam.xi.oc_sc(2) = 0.0053901;
beam.xi.oc_sc(3) = 0.0040261;

% config with sc_oc (piezo 2)
beam.nf.sc_oc(1) = 20.23; %[Hz]
beam.nf.sc_oc(2) = 117.7; %[Hz]
beam.nf.sc_oc(3) = 314.32; %[Hz]

% config with oc_oc (both piezo)
beam.nf.oc_oc(1) = 20.803; %[Hz]
beam.nf.oc_oc(2) = 117.75; %[Hz]

beam.xi.sc_oc(1) = 0.0094395;
beam.xi.sc_oc(2) = 0.0038633;
beam.xi.sc_oc(3) = 0.0038445;

beam.k.k1 = sqrt ((beam.nf.oc_sc .^2 - beam.nf.sc_sc.^2) ./ beam.nf.sc_sc .^2);
beam.k.k2 = sqrt ((beam.nf.sc_oc .^2 - beam.nf.sc_sc.^2) ./ beam.nf.sc_sc .^2);

beam.Cp.C10 = 8.3e-8;
beam.Cp.C11 = beam.Cp.C10 ./ (1 + beam.k.k1(1).^2);
beam.Cp.C12 = beam.Cp.C11 ./ (1 + beam.k.k2(1).^2);

beam.Cp.C22 = 7.46e-8;
beam.Cp.C21 = beam.Cp.C22 .* (1 + beam.k.k2(2).^2);

