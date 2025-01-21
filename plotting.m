
% plot of the modes
figure(1)
plot(x/L,phi(:,1),'k-','linewidth',2)
hold on
plot(x/L,phi(:,2),'k--','linewidth',2)
plot(x/L,phi(:,3),'k-.','linewidth',2)
plot(x/L,phi(:,4),'k:','linewidth',2)
plot(x/L,phi(:,5),'k-','linewidth',2)
xlabel('x/L')
ylabel('Modal Amplitude')
legend('1st modeshape','2nd modeshape','3rd modeshape','4th modeshape','5th modeshape')
grid

% plot of the FRF
figure;
subplot(2,1,1)
semilogy(omegav/2/pi,abs(Hsci),'--')
hold on
semilogy(omegav/2/pi,abs(Hoci),'--')
semilogy(omegav/2/pi,abs(Hsc),'--','linewidth',1)
semilogy(omegav/2/pi,abs(Hoc),'linewidth',1)
xlabel('Frequency [Hz]')
ylabel('|H| [m/N]')
grid on
axis tight
subplot(2,1,2)
plot(omegav/2/pi,angle(Hsci))
hold on
plot(omegav/2/pi,angle(Hoci),'--')
plot(omegav/2/pi,angle(Hsc),'linewidth',1)
plot(omegav/2/pi,angle(Hoc),'--','linewidth',1)
grid on
axis tight

% figure
% tiledlayout(2,1)
% ax1 = nexttile;
% semilogy(omegav/2/pi,abs(Hsc))
% hold on
% semilogy(omegav/2/pi,abs(Hoc),'--')
% xlabel('Frequency [Hz]')
% ylabel('|H| [m/N]')
% grid on
% axis tight
% legend('SC','OC')
% ax2 = nexttile;
% plot(omegav/2/pi,angle(Hsc))
% hold on
% plot(omegav/2/pi,angle(Hoc),'--')
% grid on
% axis tight
% linkaxes([ax1 ax2],'x')
% xlabel('Frequency [Hz]')
% ylabel('\phi [rad]')