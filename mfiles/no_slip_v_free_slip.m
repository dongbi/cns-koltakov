% Compares energy timeseries for no slip vs. free slip case
clear all; clc; close all;

% load /home/barthur/Desktop/noslip.mat
% noslip = subplot(2,1,1);
% hold all;
% Eb = Eb - Eb(1);
% Ep = Ep - Ep(1);
% Ea = Ep - Eb;
% Et = Ep + Ek;
% t = t/10;
% plot(t,Ep,'k--',t,Eb,'b',t,Ea,'r',t,Ek,'k:',t,Et,'k');
% plot([min(t) max(t)],[0 0],'k');
% set(gca,'XTickLabel',[]);
% ylabel('E-E_0');
% title('No-slip');
% legend('E_p','E_b','E_a','E_k','E_T','Location','Northwest');
% axis([min(t) max(t) -.7 1]);

% load /home/barthur/Desktop/freeslip.mat
% t = t/10;
% freeslip = subplot(2,1,2);
% hold all;
% Eb = Eb - Eb(1);
% Ep = Ep - Ep(1);
% Ea = Ep - Eb;
% Et = Ep + Ek;
% plot(t,Ep,'k--',t,Eb,'b',t,Ea,'r',t,Ek,'k:',t,Et,'k');
% plot([min(t) max(t)],[0 0],'k');
% xlabel('t/T');
% ylabel('E-E_0');
% title('Free-slip');
% axis([min(t) max(t) -.7 1]);

% load /home/barthur/Desktop/freeslip.mat
% Et = Ep+Ek;
% tf = t;
% dEtdt = (Et(3:end)-Et(1:end-2))/0.01;
% tf = tf/10;
% 
% load /home/barthur/Desktop/eps2.mat
% epsilon_noslip = epsilon_tot;
%  
% load /home/barthur/Desktop/eps3.mat
% epsilon_freeslip = epsilon_tot;
% 
% figure;
% hold all;
% t = t/10;
% plot(tf(2:end-1),dEtdt,'k',t,epsilon_noslip,'b',t,epsilon_freeslip,'r');
% xlabel('t/T');
% ylabel('Energy flux');
% title('No-slip');
% legend('dEt/dt','\epsilon no slip','\epsilon free slip')

load /home/barthur/Desktop/noslip.mat
Ebno = Eb - Eb(1);
Epno = Ep - Ep(1);
Eano = Ep - Eb;
Eano = Eano - Eano(1);
Ekno = Ek;
Etno = Ep + Ek;
load /home/barthur/Desktop/freeslip.mat
Ebfree = Eb - Eb(1);
Epfree = Ep - Ep(1);
Eafree = Ep - Eb;
Eafree = Eafree - Eafree(1);
Ekfree = Ek;
Etfree = Ep + Ek;

figure;
subplot(4,1,1)
plot(t,Epfree-Epno,'k');
subplot(4,1,2)
plot(t,Ebfree-Ebno,'k');
subplot(4,1,3)
plot(t,Eafree-Eano,'k');
subplot(4,1,4)
plot(t,Ekfree-Ekno,'k');
xlabel('t [s]');

figure;
subplot(4,1,1)
plot(t,Epno,'k',t,Epfree,'k--');
axis tight;
ylabel('E_p-E_p_0');
legend('no slip','free slip','Location','NorthWest');
subplot(4,1,2)
plot(t,Ebno,'k',t,Ebfree,'k--');
axis tight;
ylabel('E_b-E_b_0');
subplot(4,1,3)
plot(t,Eano,'k',t,Eafree,'k--');
axis tight;
ylabel('E_a-E_a_0');
subplot(4,1,4)
plot(t,Ekno,'k',t,Ekfree,'k--');
axis tight;
ylabel('E_k-E_k_0');
xlabel('t [s]');


