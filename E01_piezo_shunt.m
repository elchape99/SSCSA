clear all
close all
clc

% Dati trave
beam = load('Data.mat');
nf=beam.nf1_sc;
csit=beam.csi1;

% load('beam.mat')
% nf=[22 138 385 754];
% csit=[0.38 0.28 0.32 0.35]./100;

% nf=[22 28 385 754];
% csit=[0.45 0.28 0.32 0.35]./100;

% nf=[22 23 385 754];
% csit=[0.45 0.38 0.32 0.35]./100;


nmodes=length(nf);
omi=2*pi*nf; % natural pulse of the modes, no action of the piezo 

k=[0.27 0.1 0.07 0.12];
omoct=omi.*sqrt(1+k.^2); % natural freq open circuit

L = 0.30;   % length [m]
b = 0.015;  % width  [m]
h = 0.002; % height [m]

le = 0.005; 
ne = round(L/le);
a=le/2* ones(1,2*ne);
ltot=2*a*ne;

% discretize the lenght of the beam
x = (0:ltot(1)/ne:ltot(1)*(1-1/ne));
[erre,ci]=size(beam.PHI);
phi=zeros(60,6);

for ii=1:ci
    % discretize the magnitudo of mode shape PHI 
    phi(:,ii)=beam.PHI(:,ii)./max(abs(beam.PHI(:,ii)));
end

% find the Frequency Respone Function (FRF)
m=length(beam.PHI(:,1)); % measuring point, in our case equal to focing point
omegav=2*pi*[0:0.01:nf(end)*1.3].'; % system excitation frequency !!!

[Hsci, Hsc] = ComputeFRF (beam.PHI(:,1:4), beam.nf1_sc, omegav, beam.csi1, m, m);
[Hoci, Hoc] = ComputeFRF (beam.PHI(:,1:4), beam.nf1_oc, omegav, beam.csi1, m, m);

% plot all the staff
plotting

% tuning R and LR
c0 = 33.26*10^-9; % measured capacity !!!!!
% equivalent capacitance related at each mode(cpi into slides)
cpt=[c0/(1+k(1)^2) c0/(1+k(1)^2)/(1+k(2)^2) c0/(1+k(1)^2)/(1+k(2)^2)/(1+k(3)^2) c0/(1+k(1)^2)/(1+k(2)^2)/(1+k(3)^2)/(1+k(4)^2)];
rt=zeros(nmodes,1);

for jj = 1:nmodes
    Hr=zeros(length(omegav),1);
    Hlr=zeros(length(omegav),1);
    omf=sqrt((omoct(jj)^2+omi(jj)^2)/2);
    ta=1/omf;
    rt(jj)=ta/cpt(jj);
    omet(jj)=omoct(jj);
    lt(jj)=1/(omet(jj)^2*cpt(jj));
    xiet(jj)=sqrt(3)/2*sqrt((omoct(jj)^2-omi(jj)^2)/(omoct(jj)^2+omi(jj)^2));
    rlt(jj)=2*xiet(jj)*sqrt(lt(jj)/cpt(jj));
    for ii=1:nmodes
        r=rt(jj);
        l=lt(jj);
        cpi=cpt(ii);
        tau=r*cpi;
        ome=sqrt(1/(l*cpi));
        xie=rlt(jj)/2*sqrt(cpi/l);
        csi=csit(ii);
        omz=omi(ii);
        omoc=omoct(ii);
        fi_m=beam.PHI(m,ii);
        fi_f=beam.PHI(m,ii);
        % !!!! chiedere questa formul aalla prof perche è diversa dalla
        % teoria , in georia non ce fi_m, fi_f
        Hri=fi_m*fi_f*((1+1i*omegav*tau)./(-omegav.^2.*(1+2*csi*omz*tau)+1i*omegav.*(2*csi.*omz+omoc.^2*tau-tau*omegav.^2)+omz^2)); 
        Hrli=fi_m*fi_f*(-omegav.^2+ome^2+2*1i*xie*ome*omegav)./(omegav.^4-omegav.^2.*(ome^2+4*csi*xie*omz*ome+omoc^2)+omz^2*ome^2+1i*omegav.*(2*xie*ome*(omoc^2-omegav.^2)+2*csi*omz*(ome^2-omegav.^2)));
        Hr=Hri+Hr;
        Hlr=Hrli+Hlr;
    end
    figure
    ax1=subplot(2,1,1)
    semilogy(omegav/2/pi,abs(Hsc),'--')
    hold on
    semilogy(omegav/2/pi,abs(Hoc),'--')
    semilogy(omegav/2/pi,abs(Hr))
    semilogy(omegav/2/pi,abs(Hlr))


    xlabel('Frequency [Hz]')
    ylabel('|H| [m/N]')
    grid on
    axis tight
    legend('SC','OC','R','RL')

    ax2=subplot(2,1,2)
    plot(omegav/2/pi,angle(Hsc),'--')
    hold on
    plot(omegav/2/pi,angle(Hoc),'--')
    plot(omegav/2/pi,angle(Hr))
    plot(omegav/2/pi,angle(Hlr))

    grid on
    axis tight
    xlabel('Frequency [Hz]')
    ylabel('\phi [rad]')

    linkaxes([ax1,ax2],'x');


    % figure
    % subplot(2,1,1)
    % plot(omegav/2/pi,abs(Hsc),'--')
    % hold on
    % plot(omegav/2/pi,abs(Hoc),'--')
    % plot(omegav/2/pi,abs(Hr))
    % plot(omegav/2/pi,abs(Hlr))
    % xlabel('Frequency [Hz]')
    % ylabel('|H| [m/N]')
    % grid on
    % axis tight
    % legend('SC','OC','R','RL')
    % subplot(2,1,2)
    % plot(omegav/2/pi,angle(Hsc),'--')
    % hold on
    % plot(omegav/2/pi,angle(Hoc),'--')
    % plot(omegav/2/pi,angle(Hr))
    % plot(omegav/2/pi,angle(Hlr))
    % grid on
    % axis tight
    % xlabel('Frequency [Hz]')
    % ylabel('\phi [rad]')
    % clear Hr Hlr
end

