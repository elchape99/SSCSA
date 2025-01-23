clear all
close all

[file,path] = uigetfile('MultiSelect','on');
cd(path);

for ii=1:length(file)
    if iscell(file)==0
    nome=file;
    else
    nome=str2mat(file(ii));
    end
load(nome);
Freq_exp=dev7645.imps.sample{1, 1}.frequency;
Cp_exp=dev7645.imps.sample{1, 1}.param1;
R_exp=dev7645.imps.sample{1, 1}.param0;
%%%aggiunta stefano
% CnF=(CnF+0.55)*10^-9;
% cnnnn=[CnF(1:8); CnF(10:31); CnF(33:39); CnF(41:end)];
% CnF=cnnnn;
% ffffff=[fHz(1:8) ;fHz(10:31); fHz(33:39) ;fHz(41:end)];
% fHz=ffffff;
%%%%


figure(1)
hold on
plot(Freq_exp,Cp_exp*10^9)
grid on
title('Experimental measurements')
xlabel('Frequency [Hz]')
ylabel('Cp [nF]')
axis tight
% set(gca, 'XScale', 'log');

figure(2)
hold on
plot(Freq_exp,R_exp/10^3)
grid on
title('Experimental measurements')
xlabel('Frequency [Hz]')
ylabel('Resistance [k\Omega]')
axis tight
set(gca, 'XScale','log','YScale','log')
end