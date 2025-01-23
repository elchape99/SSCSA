% bar definition
% - first piezo is the one on the bottom
% - second piezo id the one on the top
clc
clear
close all

% config with sc_sc
beam.nf.sc_sc(1) = 20.122; %[Hz]
beam.nf.sc_sc(2) = 115.74; %[Hz]
beam.nf.sc_sc(3) = 314.24; %[Hz]

beam.xi.sc_sc(1) = 0.010033;
beam.xi.sc_sc(2) = 0.0047168;
beam.xi.sc_sc(3) = 0.0028293;

% config with oc_sc
beam.nf.oc_sc(1) = 20.723; %[Hz]
beam.nf.oc_sc(2) = 115.94; %[Hz]
beam.nf.oc_sc(3) = 314.81; %[Hz]

beam.xi.oc_sc(1) = 0.0086232;
beam.xi.oc_sc(2) = 0.0053901;
beam.xi.oc_sc(3) = 0.0040261;

% config with sc_oc
beam.nf.sc_oc(1) = 20.23; %[Hz]
beam.nf.sc_oc(2) = 117.7; %[Hz]
beam.nf.sc_oc(3) = 314.32; %[Hz]

beam.xi.sc_oc(1) = 0.0094395;
beam.xi.sc_oc(2) = 0.0038633;
beam.xi.sc_oc(3) = 0.0038445;
