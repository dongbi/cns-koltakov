%To analyze energy output for solitary wave cases for GRS/GRC

load /home/barthur/zang/GRC/diss-no-1.mat
% load ./energy1.mat
%Determine breaking point
[~,Ipks] = findpeaks(Ea);
[~,Ipks2] = findpeaks(-Ea);
Ib = Ipks(1);
Iend = Ipks2(2);

%Plot time series
figure(1);
hold on;
plot(t,Ep,'k--',t,Eb,'b',t,Ea,'r',t,Ek,'k-.',t,Et,'k');
plot([t(Ib) t(Ib)],[Ea(Ib) Eb(Ib)],'k:');
plot([t(Iend) t(Iend)],[Ea(Iend) Eb(Iend)],'k:');
xlabel('t [s]');
ylabel('E-E_0');
title('Time-series of Potential and Kinetic Energy');
legend('E_p','E_b','E_a','E_k','E_T','Location','EastOutside');

%Calculate bulk mixing efficiency
delta_Eb = Eb(Iend)-Eb(Ib);
delta_Ea = Ea(Iend)-Ea(Ib);
eta_b = delta_Eb/(delta_Eb-delta_Ea);
display(eta_b);

figure(2);
hold on;
load /home/barthur/zang/GRC/energy-no-1.mat
plot(t,Et,'b');
load /home/barthur/zang/GRC/energy-no-2.mat
plot(t,Et,'r');
load /home/barthur/zang/GRC/energy-no-3.mat
plot(t,Et,'k');
load /home/barthur/zang/GRC/energy-no-4.mat
plot(t,Et,'m');
load /home/barthur/zang/GRC/energy-free-1.mat
plot(t,Et,'b--');
load /home/barthur/zang/GRC/energy-free-2.mat
plot(t,Et,'r--');
load /home/barthur/zang/GRC/energy-free-3.mat
plot(t,Et,'k--');
load /home/barthur/zang/GRC/energy-free-4.mat
plot(t,Et,'m--');
xlabel(E_b-E_b_0);
ylabel('t [s]');
hold off;

% figure(3);
% hold on;
% load /home/barthur/zang/GRC/energy-no-3.mat
% plot(t,Eb,'b',t,Ea,'r');
% load /home/barthur/zang/GRC/energy-free-3.mat
% plot(t,Eb,'b--',t,Ea,'r--');
% xlabel(E-E_0);
% ylabel('t [s]');
% hold off;

figure(4);
hold on;
load /home/barthur/zang/GRC/energy-no-1.mat
%Determine breaking point
[~,Ipks] = findpeaks(Ea);
Ib = Ipks(1);
%Calculate bulk mixing efficiency
delta_Eb = Eb(Ib:end)-Eb(Ib);
delta_Ea = Ea(Ib:end)-Ea(Ib);
eta_b = delta_Eb./(delta_Eb-delta_Ea);
plot(t(Ib:end),eta_b,'b');
load /home/barthur/zang/GRC/energy-no-2.mat
%Determine breaking point
[~,Ipks] = findpeaks(Ea);
Ib = Ipks(1);
%Calculate bulk mixing efficiency
delta_Eb = Eb(Ib:end)-Eb(Ib);
delta_Ea = Ea(Ib:end)-Ea(Ib);
eta_b = delta_Eb./(delta_Eb-delta_Ea);
plot(t(Ib:end),eta_b,'r');
load /home/barthur/zang/GRC/energy-no-3
[~,Ipks] = findpeaks(Ea);
Ib = Ipks(1);
%Calculate bulk mixing efficiency
delta_Eb = Eb(Ib:end)-Eb(Ib);
delta_Ea = Ea(Ib:end)-Ea(Ib);
eta_b = delta_Eb./(delta_Eb-delta_Ea);
plot(t(Ib:end),eta_b,'k');
load /home/barthur/zang/GRC/energy-no-4.mat
%Determine breaking point
[~,Ipks] = findpeaks(Ea);
Ib = Ipks(1);
%Calculate bulk mixing efficiency
delta_Eb = Eb(Ib:end)-Eb(Ib);
delta_Ea = Ea(Ib:end)-Ea(Ib);
eta_b = delta_Eb./(delta_Eb-delta_Ea);
plot(t(Ib:end),eta_b,'m');
load /home/barthur/zang/GRC/energy-free-1.mat
%Determine breaking point
[~,Ipks] = findpeaks(Ea);
Ib = Ipks(1);
%Calculate bulk mixing efficiency
delta_Eb = Eb(Ib:end)-Eb(Ib);
delta_Ea = Ea(Ib:end)-Ea(Ib);
eta_b = delta_Eb./(delta_Eb-delta_Ea);
plot(t(Ib:end),eta_b,'b--');
load /home/barthur/zang/GRC/energy-free-2.mat
%Determine breaking point
[~,Ipks] = findpeaks(Ea);
Ib = Ipks(1);
%Calculate bulk mixing efficiency
delta_Eb = Eb(Ib:end)-Eb(Ib);
delta_Ea = Ea(Ib:end)-Ea(Ib);
eta_b = delta_Eb./(delta_Eb-delta_Ea);
plot(t(Ib:end),eta_b,'r--');
load /home/barthur/zang/GRC/energy-free-3
[~,Ipks] = findpeaks(Ea);
Ib = Ipks(1);
%Calculate bulk mixing efficiency
delta_Eb = Eb(Ib:end)-Eb(Ib);
delta_Ea = Ea(Ib:end)-Ea(Ib);
eta_b = delta_Eb./(delta_Eb-delta_Ea);
plot(t(Ib:end),eta_b,'k--');
load /home/barthur/zang/GRC/energy-free-4.mat
%Determine breaking point
[~,Ipks] = findpeaks(Ea);
Ib = Ipks(1);
%Calculate bulk mixing efficiency
delta_Eb = Eb(Ib:end)-Eb(Ib);
delta_Ea = Ea(Ib:end)-Ea(Ib);
eta_b = delta_Eb./(delta_Eb-delta_Ea);
plot(t(Ib:end),eta_b,'m--');
axis([0 40 0 1]);
hold off;

%%
figure(5);
hold on;
load /home/barthur/zang/GRC/energy-no-1.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt,'b','LineWidth',2);
load /home/barthur/zang/GRC/energy-no-2.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt,'Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/energy-no-3.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt,'r','LineWidth',2);
% load /home/barthur/zang/GRC/energy-no-4.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% dEadt = (Ea(3:end)-Ea(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt,'m');
load /home/barthur/zang/GRC/energy-free-1.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt,'b--','LineWidth',2);
load /home/barthur/zang/GRC/energy-free-2.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt,'--','Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/energy-free-3.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt,'r--','LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-4.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% dEadt = (Ea(3:end)-Ea(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt,'m--');
axis([0 40 0 1e-5]);
xlabel('t [s]','FontSize',16);
ylabel('dE_b/dt [m^2/s^3]','FontSize',16);
legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',16,'DataAspectRatio',[1 5*10^-7 1]);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
box on;
hold off;

%%
figure(6);
hold on;
load /home/barthur/zang/GRC/diss-no-1.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'b','LineWidth',2);
load /home/barthur/zang/GRC/diss-no-2.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/diss-no-3.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'r','LineWidth',2);
% load /home/barthur/zang/GRC/energy-no-4.mat
% plot(t,epsilon_tot,'m');
load /home/barthur/zang/GRC/diss-free-1.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'b--','LineWidth',2);
load /home/barthur/zang/GRC/diss-free-2.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'--','Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/diss-free-3.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'r--','LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-4.mat
% plot(t,epsilon_tot,'m--');
% axis([0 40 0 1e-5]);
xlabel('t [s]','FontSize',16);
ylabel('\epsilon [m^2/s^3]','FontSize',16);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',16,'DataAspectRatio',[1 5*10^-7 1]);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
box on;
hold off;

%%
figure(7);
a = 0.6;
subplot(4,1,1)
hold on;
load /home/barthur/zang/GRC/diss-no-1.mat
plot(t,Ep,':','Color',[a a a],'LineWidth',2);
plot(t,Eb,'k','LineWidth',2);
plot(t,Ea,'--','Color',[a a a],'LineWidth',2)
plot(t(1:end-1),Ek(1:end-1),'-.','Color',[a a a],'LineWidth',2)
plot(t(1:end-1),Et(1:end-1),'Color',[a a a],'LineWidth',2);
ylabel('E-E_0 [m^2/s^2]','FontSize',14);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',14,'Xticklabel',[]);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
axis([0 40 -.5e-4 1.1e-4]);
box on;
hold off;

subplot(4,1,2)
hold on;
load /home/barthur/zang/GRC/energy-no-1.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt,'b','LineWidth',2);
load /home/barthur/zang/GRC/energy-no-2.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt,'Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/energy-no-3.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt,'r','LineWidth',2);
% load /home/barthur/zang/GRC/energy-no-4.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% dEadt = (Ea(3:end)-Ea(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt,'m');
load /home/barthur/zang/GRC/energy-free-1.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt,'b--','LineWidth',2);
load /home/barthur/zang/GRC/energy-free-2.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt,'--','Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/energy-free-3.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt,'r--','LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-4.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% dEadt = (Ea(3:end)-Ea(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt,'m--');
axis([0 40 0 1e-5]);
ylabel('dE_b/dt [m^2/s^3]','FontSize',14);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',14,'Xticklabel',[]);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
axis([0 40 1.1e-6 10e-6]);
box on;
hold off;

subplot(4,1,3)
hold on;
load /home/barthur/zang/GRC/diss-no-1.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'b','LineWidth',2);
load /home/barthur/zang/GRC/diss-no-2.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/diss-no-3.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'r','LineWidth',2);
% load /home/barthur/zang/GRC/energy-no-4.mat
% plot(t,epsilon_tot,'m');
load /home/barthur/zang/GRC/diss-free-1.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'b--','LineWidth',2);
load /home/barthur/zang/GRC/diss-free-2.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'--','Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/diss-free-3.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'r--','LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-4.mat
% plot(t,epsilon_tot,'m--');
% axis([0 40 0 1e-5]);
ylabel('\epsilon [m^2/s^3]','FontSize',14);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',14,'Xticklabel',[]);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
axis([0 40 0 7e-6]);
box on;
hold off;

subplot(4,1,4)
hold on;
load /home/barthur/zang/GRC/diss-no-1.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt./(dEbdt+epsilon_tot(2:end-1)),'b','LineWidth',2);
load /home/barthur/zang/GRC/diss-no-2.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt./(dEbdt+epsilon_tot(2:end-1)),'Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/diss-no-3.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt./(dEbdt+epsilon_tot(2:end-1)),'r','LineWidth',2);
% load /home/barthur/zang/GRC/energy-no-4.mat
% plot(t,epsilon_tot,'m');
load /home/barthur/zang/GRC/diss-free-1.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt./(dEbdt+epsilon_tot(2:end-1)),'b--','LineWidth',2);
load /home/barthur/zang/GRC/diss-free-2.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt./(dEbdt+epsilon_tot(2:end-1)),'--','Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/diss-free-3.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
plot(t_dt,dEbdt./(dEbdt+epsilon_tot(2:end-1)),'r--','LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-4.mat
% plot(t,epsilon_tot,'m--');
% axis([0 40 0 1e-5]);
xlabel('t [s]','FontSize',14);
ylabel('\eta','FontSize',14);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',14);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
axis([0 40 0.3 1]);
box on;
hold off;

%%
load /home/barthur/zang/GRC/no-noforcing.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
figure;
plot(t_dt,dEbdt,'k','LineWidth',2);
figure;
plot(t(1:end-1),epsilon_tot(1:end-1),'k','LineWidth',2);
figure;
plot(t(1:end-1),phi_d(1:end-1),'k','LineWidth',2);

%%
load /home/barthur/zang/GRC/no-noforcing.mat
dEbdt_no = (Eb(3:320)-Eb(1:318))/(2*25*0.005);
epsilon_tot_no = epsilon_tot(1:320);

figure(8);
a = 0.6;
subplot(4,1,1)
hold on;
load /home/barthur/zang/GRC/diss-no-1.mat
plot(t,Ep,':','Color',[a a a],'LineWidth',2);
plot(t,Eb,'k','LineWidth',2);
plot(t,Ea,'--','Color',[a a a],'LineWidth',2)
plot(t(1:end-1),Ek(1:end-1),'-.','Color',[a a a],'LineWidth',2)
plot(t(1:end-1),Et(1:end-1),'Color',[a a a],'LineWidth',2);
ylabel('E-E_0 [m^2/s^2]','FontSize',14);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',14,'Xticklabel',[]);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
axis([0 40 -.5e-4 1.1e-4]);
box on;
hold off;

subplot(4,1,2)
hold on;
load /home/barthur/zang/GRC/energy-no-1.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005)-dEbdt_no;
dEbdt(dEbdt<0) = 0;
plot(t_dt,dEbdt,'b','LineWidth',2);
load /home/barthur/zang/GRC/energy-no-2.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005)-dEbdt_no;
dEbdt(dEbdt<0) = 0;
plot(t_dt,dEbdt,'Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/energy-no-3.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005)-dEbdt_no;
dEbdt(dEbdt<0) = 0;
plot(t_dt,dEbdt,'r','LineWidth',2);
% load /home/barthur/zang/GRC/energy-no-4.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% dEadt = (Ea(3:end)-Ea(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt,'m');
% load /home/barthur/zang/GRC/energy-free-1.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt,'b--','LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-2.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt,'--','Color',[0 .7 0],'LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-3.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt,'r--','LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-4.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% dEadt = (Ea(3:end)-Ea(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt,'m--');
ylabel('dE_b/dt [m^2/s^3]','FontSize',14);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',14,'Xticklabel',[]);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
axis([0 40 0 4e-6]);
box on;
hold off;

subplot(4,1,3)
hold on;
load /home/barthur/zang/GRC/diss-no-1.mat
plot(t(1:end-1),epsilon_tot(1:end-1)-epsilon_tot_no(1:end-1),'b','LineWidth',2);
load /home/barthur/zang/GRC/diss-no-2.mat
plot(t(1:end-1),epsilon_tot(1:end-1)-epsilon_tot_no(1:end-1),'Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/diss-no-3.mat
plot(t(1:end-1),epsilon_tot(1:end-1)-epsilon_tot_no(1:end-1),'r','LineWidth',2);
% load /home/barthur/zang/GRC/energy-no-4.mat
% plot(t,epsilon_tot,'m');
% load /home/barthur/zang/GRC/diss-free-1.mat
% plot(t(1:end-1),epsilon_tot(1:end-1),'b--','LineWidth',2);
% load /home/barthur/zang/GRC/diss-free-2.mat
% plot(t(1:end-1),epsilon_tot(1:end-1),'--','Color',[0 .7 0],'LineWidth',2);
% load /home/barthur/zang/GRC/diss-free-3.mat
% plot(t(1:end-1),epsilon_tot(1:end-1),'r--','LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-4.mat
% plot(t,epsilon_tot,'m--');
% axis([0 40 0 1e-5]);
ylabel('\epsilon [m^2/s^3]','FontSize',14);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',14,'Xticklabel',[]);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
axis([0 40 0 7e-6]);
box on;
hold off;

subplot(4,1,4)
hold on;
load /home/barthur/zang/GRC/diss-no-1.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005)-dEbdt_no;
dEbdt(dEbdt<0) = 0;
epsilon_tot = epsilon_tot - epsilon_tot_no;
plot(t_dt,dEbdt./(dEbdt+epsilon_tot(2:end-1)),'b','LineWidth',2);
load /home/barthur/zang/GRC/diss-no-2.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005)-dEbdt_no;
dEbdt(dEbdt<0) = 0;
epsilon_tot = epsilon_tot - epsilon_tot_no;
plot(t_dt,dEbdt./(dEbdt+epsilon_tot(2:end-1)),'Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/diss-no-3.mat
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005)-dEbdt_no;
dEbdt(dEbdt<0) = 0;
epsilon_tot = epsilon_tot - epsilon_tot_no;
plot(t_dt,dEbdt./(dEbdt+epsilon_tot(2:end-1)),'r','LineWidth',2);
% load /home/barthur/zang/GRC/energy-no-4.mat
% plot(t,epsilon_tot,'m');
% load /home/barthur/zang/GRC/diss-free-1.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt./(dEbdt+epsilon_tot(2:end-1)),'b--','LineWidth',2);
% load /home/barthur/zang/GRC/diss-free-2.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt./(dEbdt+epsilon_tot(2:end-1)),'--','Color',[0 .7 0],'LineWidth',2);
% load /home/barthur/zang/GRC/diss-free-3.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt./(dEbdt+epsilon_tot(2:end-1)),'r--','LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-4.mat
% plot(t,epsilon_tot,'m--');
% axis([0 40 0 1e-5]);
xlabel('t [s]','FontSize',14);
ylabel('\eta','FontSize',14);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',14);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
axis([0 40 0 1]);
box on;
hold off;

%%
figure(9);
a = 0.6;
subplot(4,1,1)
hold on;
load /home/barthur/zang/GRC/diss-no-1.mat
plot(t,Ep,':','Color',[a a a],'LineWidth',2);
plot(t,Eb,'k','LineWidth',2);
plot(t,Ea,'--','Color',[a a a],'LineWidth',2)
plot(t(1:end-1),Ek(1:end-1),'-.','Color',[a a a],'LineWidth',2)
plot(t(1:end-1),Et(1:end-1),'Color',[a a a],'LineWidth',2);
ylabel('E-E_0 [m^2/s^2]','FontSize',14);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',14,'Xticklabel',[]);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
axis([0 40 -.5e-4 1.1e-4]);
box on;
hold off;

subplot(4,1,2)
hold on;
load /home/barthur/zang/GRC/phid-no-1.mat
phi_d = phi_d(1:319);
plot(t(1:319),phi_d,'b','LineWidth',2);
load /home/barthur/zang/GRC/phid-no-2.mat
phi_d = phi_d(1:319);
plot(t(1:319),phi_d,'Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/phid-no-3.mat
phi_d = phi_d(1:319);
plot(t(1:319),phi_d,'r','LineWidth',2);
% load /home/barthur/zang/GRC/energy-no-4.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% dEadt = (Ea(3:end)-Ea(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt,'m');
load /home/barthur/zang/GRC/phid-free-1.mat
phi_d = phi_d(1:319);
plot(t(1:319),phi_d,'b--','LineWidth',2);
load /home/barthur/zang/GRC/phid-free-2.mat
phi_d = phi_d(1:319);
plot(t(1:319),phi_d,'--','Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/phid-free-3.mat
phi_d = phi_d(1:319);
plot(t(1:319),phi_d,'r--','LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-4.mat
% dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
% dEadt = (Ea(3:end)-Ea(1:end-2))/(2*25*0.005);
% plot(t_dt,dEbdt,'m--');
axis([0 40 0 1e-5]);
ylabel('dE_b/dt [m^2/s^3]','FontSize',14);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',14,'Xticklabel',[]);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
axis([0 40 1.1e-7 10e-7]);
box on;
hold off;

subplot(4,1,3)
hold on;
load /home/barthur/zang/GRC/diss-no-1.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'b','LineWidth',2);
load /home/barthur/zang/GRC/diss-no-2.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/diss-no-3.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'r','LineWidth',2);
% load /home/barthur/zang/GRC/energy-no-4.mat
% plot(t,epsilon_tot,'m');
load /home/barthur/zang/GRC/diss-free-1.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'b--','LineWidth',2);
load /home/barthur/zang/GRC/diss-free-2.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'--','Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/diss-free-3.mat
plot(t(1:end-1),epsilon_tot(1:end-1),'r--','LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-4.mat
% plot(t,epsilon_tot,'m--');
% axis([0 40 0 1e-5]);
ylabel('\epsilon [m^2/s^3]','FontSize',14);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',14,'Xticklabel',[]);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
axis([0 40 0 7e-6]);
box on;
hold off;

subplot(4,1,4)
hold on;
load /home/barthur/zang/GRC/phid-no-1.mat
phi_d = phi_d(2:319);
load /home/barthur/zang/GRC/diss-no-1.mat
plot(t_dt,phi_d./(phi_d+epsilon_tot(2:end-1)),'b','LineWidth',2);
load /home/barthur/zang/GRC/phid-no-2.mat
phi_d = phi_d(2:319);
load /home/barthur/zang/GRC/diss-no-2.mat
plot(t_dt,phi_d./(phi_d+epsilon_tot(2:end-1)),'Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/phid-no-3.mat
phi_d = phi_d(2:319);
load /home/barthur/zang/GRC/diss-no-3.mat
plot(t_dt,phi_d./(phi_d+epsilon_tot(2:end-1)),'r','LineWidth',2);
% load /home/barthur/zang/GRC/energy-no-4.mat
% plot(t,epsilon_tot,'m');
load /home/barthur/zang/GRC/phid-free-1.mat
phi_d = phi_d(2:319);
load /home/barthur/zang/GRC/diss-free-1.mat
plot(t_dt,phi_d./(phi_d+epsilon_tot(2:end-1)),'b--','LineWidth',2);
load /home/barthur/zang/GRC/phid-free-2.mat
phi_d = phi_d(2:319);
load /home/barthur/zang/GRC/diss-free-2.mat
plot(t_dt,phi_d./(phi_d+epsilon_tot(2:end-1)),'--','Color',[0 .7 0],'LineWidth',2);
load /home/barthur/zang/GRC/phid-free-3.mat
phi_d = phi_d(2:319);
load /home/barthur/zang/GRC/diss-free-3.mat
plot(t_dt,phi_d./(phi_d+epsilon_tot(2:end-1)),'r--','LineWidth',2);
% load /home/barthur/zang/GRC/energy-free-4.mat
% plot(t,epsilon_tot,'m--');
% axis([0 40 0 1e-5]);
xlabel('t [s]','FontSize',14);
ylabel('\eta','FontSize',14);
% legend('\xi=1.1','\xi=0.8','\xi=0.6');
set(gca,'FontSize',14);
% set(gca,'Ytick',10^-6*[0 1 2 3 4 5 6 7 8 9 10]);
% set(gca,'YtickLabel',[0 1 2 3 4 5 6 7 8 9 10]);
% text(-4,1e-5,'(\times10^{-6})','FontSize',14);
axis([0 40 0 1]);
box on;
hold off;