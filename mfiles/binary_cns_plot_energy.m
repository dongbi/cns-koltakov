% Plots CNS simulation data from .mat file created with
% binary_cns_load.m
clear all; clc; close all;

load energy1.mat
% load /home/barthur/Desktop/energy.mat
% load /home/barthur/Desktop/solitary_energy.mat

figure(1);
plot(t,Ep,'k--',t,Eb,'b',t,Ea,'r',t,Ek,'k:',t,Et,'k');
xlabel('t [s]');
ylabel('E-E_0');
title('Time-series of Potential and Kinetic Energy');
legend('E_p','E_b','E_a','E_k','E_T','Location','EastOutside');

figure(2)
plot(t_dt,dEpdt,'k--',t_dt,dEbdt,'b',t_dt,dEadt,'r',t_dt,dEkdt,'k:',t_dt,dEtdt,'k', ...
    t,F_Ep,'g',t,F_Ek,'m',t_load,epsilon_tot,'c',t_load,epsilon_tke,'y',t,phi_d,'b--');
xlabel('t [s]');
ylabel('Energy flux [W]');
title('Time-series of Energy Fluxes');

% figure(3)
% F_pe = interp(F_pe,100)/1000; F_pe = F_pe(2:end-1);
% F_ke = interp(F_ke,100)/1000; F_ke = F_ke(2:end-1);
% epsilon_tot = interp(epsilon_tot,100); epsilon_tot = epsilon_tot(2:end-1);
% plot(t_dt,dEpdt+dEkdt-F_pe-F_ke+epsilon_tot);

% Eb_net = zeros(1,length(t)-1);
% for n=2:length(t)
%     Eb_net(n-1) = Eb(n) + F_Eb(n-1)*0.1;
% end
% figure(4);
% plot(t(2:end),Eb_net,'k',t,Eb,'b');
% 
% 
% plot(t,F_Ep,'k',t,F_Eb,'b')
% plot(t_dt,dEpdt,'k--',t,-F_Ep,'k',t_dt,dEbdt,'b--',t,F_Eb,'b')
dEbdt = (Eb(3:end)-Eb(1:end-2))/(.1*2); 
plot(t,phi_d,'b',t,epsilon_tot,'k',t_dt,dEbdt,'r');    
    