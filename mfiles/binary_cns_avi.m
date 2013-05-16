%Prints images from .figs saved with binary_cns_isoplot_movie.m
%In terminal, create .avi (or other) with:
%$ avconv -f image2 -i ./isoplot%d.jpg -s vga ./isoplot_movie.avi

clear all; close all; clc;

addpath '/home/barthur/zang/isoplot_movie_2/';

Nframes = 296;
for n=1:Nframes;
    display(n);
    openfig(['isoplot_',num2str(n),'fig.fig']);
    print('-djpeg','-r300',['/home/barthur/zang/isoplot_movie_2/isoplot',num2str(n)]);
    close(gcf);
end

% count = 1;
% Nframes = 20;
% 
% for n=1:Nframes
%     openfig(['isoplot_',num2str(n),'fig.fig']);
%     
%     if count<=Nframes;
%         f = getframe(gcf);
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
% end
% 
% imwrite(im,map,movie_name,'DelayTime',0.1,'LoopCount',inf);


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
