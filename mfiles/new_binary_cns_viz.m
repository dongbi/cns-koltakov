% Displays CNS simulation data from binary output files
clear all; close all; clc;
directory = '/home/barthur/zang/2D_test/';

% PLOTTING OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timestep = 10; %actual timestep
delta_ts = 0; %averaging
save_timestep_period = 10;
display_grid = false;
display_density = true;
display_scalar = true;
display_pressure = false;
display_averages = false;
display_lambda2 = false;
plot_velocity_z_slice=true; %x-y plot
plot_velocity_x_slice=false; %z-y plot
plot_velocity_y_slice=false; %x-z plot
show_grid_lines = false;

% global grid dims
g_ni = 128;
g_nj = 32; 
g_nk = 1;

% # procs in each dim
npx=4;
npy=1;
npz=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local grid dims
% halo <- 1
l_ni_h1 = g_ni/npx + 2;
l_nj_h1 = g_nj/npy + 2;
l_nk_h1 = g_nk/npz + 2;
% halo <- 2
l_ni_h2 = g_ni/npx + 4;
l_nj_h2 = g_nj/npy + 4;
l_nk_h2 = g_nk/npz + 4;

total_size_h1 = l_ni_h1 * l_nj_h1 * l_nk_h1;  % flat array_3d size w/ 1 halo cell
total_size_h2 = l_ni_h2 * l_nj_h2 * l_nk_h2;  % flat array_3d size w/ 2 halo cells
int_sz = 4;
dbl_sz = 8;

z_slice = round(g_nk/2);
x_slice = round(g_ni/2);
y_slice = round(g_nj); %1;

%adjust timestep
timestep = timestep - save_timestep_period;
timestep = timestep/save_timestep_period;

% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xp,Yp,Zp] = load_binary_grid(directory,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);

if(display_grid)    
    figure;
    for j=1:g_nj
        line(Xp(:,j,z_slice),Yp(:,j,z_slice), 'Color', 'black');     
    end
    for i=1:g_ni
        line(Xp(i,:,z_slice),Yp(i,:,z_slice), 'Color', 'black');     
    end      
    axis('image');
    title('Computational Grid: longitudal $$z$$-slice','Interpreter','LaTex','FontSize',16);
end

% DENSITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_density)
   [rho] = load_binary_density(directory,timestep,...
       delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
   figure;
   pcolor(Xp,Yp,rho);
   colorbar;
   axis equal;
end

% VELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(plot_velocity_x_slice || plot_velocity_y_slice || plot_velocity_z_slice)
    [U,V,W] = load_binary_velocity(directory,timestep,...
        delta_ts/save_timestep_period, npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
     
    % PLOT VELOCITY: Z-SLICE
    if(plot_velocity_z_slice)            
        figure;    
%         xz = squeeze(Xp(:,:,z_slice)); yz = squeeze(Yp(:,:,z_slice)); 
%         uz = squeeze(U(:,:,z_slice));  vz = squeeze(V(:,:,z_slice)); wz = squeeze(W(:,:,z_slice));         
%         uz = uz - mean(mean(uz));
%         vz = vz - mean(mean(vz));
%         wz = wz - mean(mean(wz));
%         cav = abs(curl(yz,xz,vz,uz));
%         f = cav;
%         h=pcolor(xz, yz, f); set( h, 'EdgeColor', 'none');
%         %contour(xz, yz, f, 50); 
%         colorbar;
%         %colormap bone;caxis([200 max(max(f))-min(min(f))])
%         hold on;
%         quiver(xz, yz, uz, vz, 'Color', [0 0 0]);         
%         %streamslice(yz,xz, uz,vz, 10); 
%         %bathymetry    
%         line(xz(:,1),yz(:,1), 'Color', 'black', 'LineWidth', 2);
%         axis('image');
%         set(gca, 'FontSize', 14);   
%         xlabel('$$X$$','Interpreter','LaTex','FontSize',16);     
%         ylabel('$$Y$$','Interpreter','LaTex','FontSize',16);   
%         title('Fluid velocity $$(u,v)$$: longitudal $$z$$-slice','Interpreter','LaTex','FontSize',16);
        pcolor(Xp,Yp,V);
        colorbar;
        axis equal;
    end
    
    % PLOT VELOCITY: X-SLICE
    if(plot_velocity_x_slice)    
        zx = squeeze(Zp(x_slice,:,:)); yx = squeeze(Yp(x_slice,:,:)); 
        ux = squeeze(U(x_slice,:,:));  vx = squeeze(V(x_slice,:,:)); wx = squeeze(W(x_slice,:,:)); 
        ux = ux - mean(mean(ux));
        vx = vx - mean(mean(vx));
        wx = wx - mean(mean(wx));        
        figure;    
        h=pcolor(zx, yx, ux); set( h, 'EdgeColor', 'none');
        %contour(zx, yx, ux, 50); 
        colorbar;
        hold on;
        quiver(zx, yx, wx, vx, 'Color', [0 0 0]);            
        %bathymetry    
        line(zx(1,:),yx(1,:), 'Color', 'black', 'LineWidth', 2);
        axis('image');
        set(gca, 'FontSize', 14);   
        xlabel('$$Z$$','Interpreter','LaTex','FontSize',16);     
        ylabel('$$Y$$','Interpreter','LaTex','FontSize',16);   
        title('Fluid velocity $$(w,v)$$: vertical $$x$$-slice','Interpreter','LaTex','FontSize',16);        
    end  
    
    % PLOT VELOCITY: Y-SLICE
    if(plot_velocity_y_slice)    
        xy = squeeze(Xp(:,y_slice,:)); zy = squeeze(Zp(:,y_slice,:)); 
        depth = abs(squeeze(Yp(:,1,:)));
        uy = squeeze(U(:,y_slice,:));  vy = squeeze(V(:,y_slice,:)); wy = squeeze(W(:,y_slice,:)); 
        uy = uy - mean(mean(uy));
        vy = vy - mean(mean(vy));
        wy = wy - mean(mean(wy));
        figure;  
        hold on;
        %h=pcolor(xy, zy, vy); set( h, 'EdgeColor', 'none');
        %colorbar;        
        %m=max(max(vy));contour(xy, zy, vy, [.1*m:.05*m:m]);
        %h=pcolor(xy, zy, squeeze(P(:,y_slice,:))); set( h, 'EdgeColor', 'none');        
        div = divergence(zy,xy,wy,uy);
        cav = abs(curl(zy,xy,wy,uy));
        %f = depth.*div; 
        %f = cav;
        f = vy;        
        %f = uy;
        h=pcolor(xy,zy,f);  set( h, 'EdgeColor', 'none');
        colorbar;
        %caxis([0 max(max(max(f)))])
        caxis([min(min(min(f))) max(max(max(f)))])
        
        contour(xy,zy,f,[0 0], 'Color', 'black', 'LineWidth', 2);
        
        %uy = uy - mean(mean(uy));        
        %wy = wy - mean(mean(wy));
        %h=pcolor(xy,zy,VORT);  set( h, 'EdgeColor', 'none');
        %contour(xy, zy, sqrt(uy.^2+wy.^2),[0.001:.001:.02]);
        quiver(xy, zy, uy, wy, 'Color', [0 0 0]);    
        %streamslice(zy,xy,uy,wy,10)                
        axis('image');
        set(gca, 'FontSize', 14);   
        xlabel('$$X$$','Interpreter','LaTex','FontSize',16);     
        ylabel('$$Z$$','Interpreter','LaTex','FontSize',16);   
        title('Fluid velocity rms $$(u,v)$$: horizontal $$y$$-slice','Interpreter','LaTex','FontSize',16);
    end     
end

% SCALAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_scalar)
   [phi] = load_binary_scalar(directory,timestep,...
       delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
   figure;
   pcolor(Xp,Yp,phi);
   colorbar;
   axis equal;
end

% PRESSURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_pressure)    
    [P] = load_binary_pressure(directory, timestep,delta_ts, npx,npy,npz, l_ni_h1,l_nj_h1,l_nk_h1);
    figure;
    pcolor(Xp,Yp,P);
    colorbar;
    axis equal;
end

% LAMBDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_lambda2)
    % lambda_2
    [Ux,Uz]=gradient(uy, Xp(2,1,1)-Xp(1,1,1), Zp(1,1,2)-Zp(1,1,1));
    [Wx,Wz]=gradient(wy, Xp(2,1,1)-Xp(1,1,1), Zp(1,1,2)-Zp(1,1,1));
    VORT = zeros(size(Ux,1),size(Ux,2));
    for i=1:size(Ux,1)
        for k=1:size(Ux,2)
            VGT = [Ux(i,k), Uz(i,k); Wx(i,k), Wz(i,k)];
            S =   .5*(VGT +  VGT');
            Omg = .5*(VGT -  VGT');
            Q = .5*(abs(trace(Omg*Omg'))-abs(trace(S*S')));
            if(Q>0)
                VORT(i,k) = Q; %1
            end
            %S2O2 = S*S + Omg*Omg;
            %ev = eig(S2O2);
            %if(ev(1) < 0 && ev(2) < 0)
            %    VORT(i,k) = 1;
            %end
        end
    end
    figure;
    h = pcolor(xy,zy,VORT); set( h, 'EdgeColor', 'none');

    figure;
    
    [curlx, curly, curlz, cav] = curl(Yp,Xp,Zp,V,U,W);
    curlmag = sqrt(curlx.^2+curly.^2+curlz.^2);    

    Pp = P-mean(mean(mean(P)));
    
    surf(xy,-depth,zy,'FaceColor','none');%[.2 .2 .2]);
    hold on;   
    pat = patch(isosurface(Xp,Yp,Zp, Pp,.0005, -curlx));
    %isonormals(Pp,pat)
    set(pat,'FaceColor','interp','EdgeColor','none');
    daspect([1,1,1])
    view([0.0, 0.02, -0.05]); axis tight
    camlight 
    lighting gouraud

    % isosurfaces of vorticity
    figure;
    surf(xy,-depth,zy,'FaceColor','none');%[.2 .2 .2]);
    hold on;
    %surf(xy,zeros(size(xy)),zy,'FaceColor','none');
    
    pat = patch(isosurface(Xp,Yp,Zp, curly,5,curlx));
    %isonormals(curlmag,pat);
    %isocolors(curlx,pat);
    set(pat,'FaceColor','interp','EdgeColor','none');
    daspect([1,1,1])
    view([0, 0.04, -0.1]); axis tight
    %view([.5, 0.1, .5]); axis tight
    camlight 
    lighting gouraud
    axis off;

    %figure;
    %h = slice(interp3(cav,Xp,Yp,Zp),.03,-.045, .05); 
    %shading interp 
    %daspect([1 1 1]); 
    %axis tight
    %colormap hot(16)
    %camlight
    %set([h(1),h(2)],'ambientstrength',.6)
end;

set(gcf, 'PaperPositionMode', 'auto');  