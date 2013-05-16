% Plots a movie of CNS simulation data from .mat file created with
% binary_cns_load.m (lateral_average=1) and saves as animated gif
clear all; clc; close all;

load /home/barthur/Desktop/1109177_all.mat

%%
% PLOTTING OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

movie = 0;
movie_name = '/home/barthur/Desktop/movie.gif';

iskip = 1;
kskip = 1;

%Figure properties
h1=figure(1);
set(h1,'Position',[100 100 1500 1000]);
set(h1,'PaperPositionMode','auto')
axis equal;
% set(gca,'DataAspectRatio',[1 .5 1]);
% axis([1.7 2 -.4 -.2]);
axis([0 x_length -z_length 0]);
xlabel('x [m]');
ylabel('z [m]');
box on;
hold on;
set(gcf,'color','w');
% set(gca,'nextplot','replacechildren','visible','off')
count = 1;
Nframes = length(t);

% for n=112
for n=1:Nframes
    str = sprintf(['n =',num2str(n),' t = ',num2str(t(n)), ...
        ' timestep = ',num2str(t(n)/dt)]);  disp(str);
    cla;
%     pcolor(x(:,:),z(:,:),squeeze(sqrt(u(:,:,n).^2+w(:,:,n).^2))); caxis([-0.2 0.2]);
    pcolor(x(:,:),z(:,:),squeeze(rho(:,:,n))); caxis([-0.015 0.015]);
%     colorbar;
    shading flat;
%     quiver(x(1:iskip:end,1:kskip:end),z(1:iskip:end,1:kskip:end), ...
%         squeeze(u(1:iskip:end,1:kskip:end,n)),squeeze(w(1:iskip:end,1:kskip:end,n)),0.1,'k');
%     axis off;
    drawnow;
    pause;
    
    if movie;
        if count<=Nframes;
            f = getframe(gcf);
            if count==1;
                [im,map] = rgb2ind(f.cdata,256,'nodither');
                im(1,1,1,Nframes) = 0;
            else
                im(:,:,1,count) = rgb2ind(f.cdata,map,'nodither');
            end
        else
            break;
        end

    count=count+1;
    end
end

if movie
    imwrite(im,map,movie_name,'DelayTime',0.1,'LoopCount',inf);
end




        
    