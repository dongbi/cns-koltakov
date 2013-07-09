%Prints images from .figs saved with binary_cns_isoplot_movie.m
%In terminal, create .avi (or other) with:
%$ avconv -f image2 -i ./isoplot%d.jpg -s vga ./isoplot_movie.avi

clear all; close all; clc;

addpath '/home/barthur/zang/isoplot_movie_2/';

Nframes = 296;
for n=1:Nframes;
    display(n);
    openfig(['isoplot_',num2str(n),'fig.fig']);
    print('-djpeg',['/home/barthur/zang/animated_gif/isoplot',num2str(n+99)]);
    close(gcf);
end

%figure properties
% h1=figure(1);
% set(h1,'Position',[100 100 1500 1000]);
% set(h1,'PaperPositionMode','auto')
% set(gcf,'color','w');
% set(gca,'nextplot','replacechildren','visible','off')
% count = 1;
% Nframes = 294;
% 
% for n=1:Nframes
%     h1 = openfig(['isoplot_',num2str(n),'fig.fig'],'new','visible');
%     set(h1,'PaperPositionMode','auto')
%     set(gcf,'color','w');
%     if count<=Nframes;
%         f = getframe(h1);
%         if count==1;
%             [im,map] = rgb2ind(f.cdata,256,'nodither');
%             im(1,1,1,Nframes) = 0;
%         else
%             im(:,:,1,count) = rgb2ind(f.cdata,map,'nodither');
%         end
%     else
%         break;
%     end
%     
%     count=count+1;
%     close(h1);
% end
% 
% imwrite(im,map,'/home/barthur/Dropbox/Public/isoplot_movie_2.gif','DelayTime',0.1,'LoopCount',Inf);


% %# figure
% % figure, set(gcf, 'Color','white')
% % Z = peaks; surf(Z);  axis tight
% % set(gca, 'nextplot','replacechildren', 'Visible','off');
% 
% %# create AVI object
% nFrames = 20;
% vidObj = VideoWriter('isoplot.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 7;
% open(vidObj);
% 
% %# create movie
% for n=1:nFrames
%    cla;
%    str = sprintf(['n = ',num2str(n)]); disp(str);
%    F = openfig(['isoplot_',num2str(n),'fig.fig']);
% %    set(gcf,'visible','on');
%    writeVideo(vidObj, getframe(F));
%    close(F);
% end
% 
% %# save as AVI file, and open it using system video player
% close(vidObj);
