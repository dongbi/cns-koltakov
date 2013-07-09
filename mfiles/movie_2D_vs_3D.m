%Create 2D vs 3D movie
clear all; close all; clc;

% load /home/barthur/Desktop/1109177_all.mat
load /home/barthur/Desktop/1109466_midslice.mat
x3 = squeeze(x);
z3 = squeeze(z);
rho3 = rho;

load christine_2D.mat
x2 = x;
z2 = y;
rho2 = rho;

%%
%Figure properties
figure(1);
set(gcf,'color','w');
set(gcf, 'InvertHardCopy', 'off');
set(gcf,'PaperUnits','normalized');
% set(gcf,'Position',[100 100 1500 1000]);
set(gcf,'PaperPositionMode','auto');

Nframes = Nt;
for n=1:Nframes;
% for n=200
    display(n);

    %3D
    h1 = axes('Position',[0.05 .5 .9 .4]);
%     h1 = subplot(2,1,1);
    hold on;
    axis off;
    axis equal;
%     set(h1,'DataAspectRatio',[1 .5 1]);
    pcolor(x3,z3,rho3(:,:,n));
    shading interp;
    hold off;

    %2D
    h2 = axes('Position',[0.05 .1 .9 .4]);
%     h2 = subplot(2,1,2);
    hold on;
    axis off;
    axis equal;
%     set(h2,'DataAspectRatio',[1 .5 1]);
    pcolor(x2,z2,rho2(:,:,n));
    shading interp;
    hold off;
    
    print('-djpeg','-r300',['/home/barthur/zang/movie_2D_vs_3D/plot',num2str(n)]);
    close(gcf);
end