% Plots a movie of density and corresponding energy
clear all; clc; close all;

load /home/barthur/zang/GRC/midslice-free-2.mat
t50 = t;
Nt50 = Nt;
load /home/barthur/zang/GRC/diss-free-2.mat
% Eb = Eb - Eb(1);
% Ep = Ep - Ep(1);
% Ea = Ep - Eb;
% Et = Ep + Ek;
% t = t/10;

Nframes = Nt50;
for n=1:Nframes;
   str = sprintf(['n =',num2str(n),' t = ',num2str(t(n)), ...
       ' timestep = ',num2str(t(n)/dt)]);  disp(str);
   
   %density
   subplot(2,1,1);
   cla;
%    contour(x(:,:),z(:,:),squeeze(rho(:,:,n)),[1, 1],'k');
   hold on;
   pcolor(x(:,:),z(:,:),squeeze(rho(:,:,n))); caxis([0.985 1.015]);
   shading flat;
   axis off;
   drawnow;
   hold off;
   
   %energy
   subplot(2,1,2)
   hold on;
   cla;
   dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*25*0.005);
   plot(t,epsilon_tot,'k',t_dt,dEbdt,'b');
%    plot(t,Ep,'k--',t,Eb,'b',t,Ea,'r',t,Ek,'k:',t,Et,'k');
%    axis([min(t) max(t) -0.25 0.25]);
%    plot(t_load/10,epsilon_tot,'b');
%    axis([min(t) max(t) 0 1e-5]);
   axis tight;
   plot([t50(n) t50(n)],[-1 1],'k');
   xlabel('t/T');
   ylabel('E-E_0');
   drawnow;
   
   pause(.02);
end