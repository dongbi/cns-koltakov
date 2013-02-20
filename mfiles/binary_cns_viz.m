% Displays CNS simulation data in binary format
clear all;
timestep = 3500;
delta_ts = 0; %averaging
save_timestep_period = 1;
timestep = timestep/save_timestep_period;
display_grid = true;
display_pressure = false;

display_averages = false;
display_lambda2 = false;

directory = '../output/';

%directory = './cos_channel_96x48x48_a10p_d5_124Kts/';
%directory = './cos_channel_96x48x48_a20p_d7_86Kts/';
%directory = './cos_channel_96x48x48_a30p_d5_102Kts/';
%directory = './sin_channel_96x48x48_200Kts/'; %d=.4

%directory = './flat_channel_96x24x48_d237_185Kts/'; %g_nj = 24;
%directory = './flat_channel_96x48x48_d474_120Kts/'; %g_nj = 48;
%directory = './flat_channel_96x48x48_d711_130Kts/'; %g_nj = 48;

%directory = './flat_channel_96x24x48_d381_125Kts/';
%directory = '../flat_channel_96x48x48_d5_200Kts/';
%directory = './flat_channel_96x48x48_d7_110Kts/';

% global grid dims
g_ni = 96;
g_nj = 48; 
g_nk = 48;

% # procs in each dim
npx=4;
npy=2;
npz=2;

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

plot_velocity_z_slice=true; z_slice = round(g_nk/2);
plot_velocity_x_slice=true; x_slice = round(g_ni/2);
plot_velocity_y_slice=true; y_slice = round(g_nj); %1;

% LOAD GRID
[Xp, Yp, Zp] = load_binary_grid(directory, npx,npy,npz, l_ni_h2,l_nj_h2,l_nk_h2);

% DISPLAY GRID
if(display_grid)    
    figure;
    for j=1:g_nj
        line(Xp(:,j,z_slice),Yp(:,j,z_slice), 'Color', 'black');     
    end
    for i=1:g_ni
        line(Xp(i,:,z_slice),Yp(i,:,z_slice), 'Color', 'black');     
    end      
    % displays grids per cpu
    %for px = 1:npx
    %    for py = 1:npy
    %        for j=l_nj_h2*(py-1)+1:l_nj_h2*(py) %g_nj+2*h
    %            line(Xp(l_ni_h2*(px-1)+1:l_ni_h2*(px),j,z_slice),Yp(l_ni_h2*(px-1)+1:l_ni_h2*(px),j,z_slice), 'Color', [(px+py)/(npx+npy) (px+py)/(npx+npy) 0]);     
    %        end
    %        for i=l_ni_h2*(px-1)+1:l_ni_h2*(px) %g_ni+2*h
    %            line(Xp(i,l_nj_h2*(py-1)+1:l_nj_h2*(py),z_slice),Yp(i,l_nj_h2*(py-1)+1:l_nj_h2*(py),z_slice), 'Color', [(px+py)/(npx+npy) (px+py)/(npx+npy) 0]);     
    %        end  
    %    end
    %end
    axis('image');
    title('Computational Grid: longitudal $$z$$-slice','Interpreter','LaTex','FontSize',16);
end

% DISPLAY PRESSURE
if(display_pressure)    
    % load P
    filename = [directory 'pressure.'];
    for pk = 1:npz
        for pj = 1:npy
            for pi = 1:npx
                proc = (pi-1)*npz*npy + (pj-1)*npz + pk-1;
                fid = fopen([filename num2str(proc)],'r');   
                
                if(fseek(fid,timestep*(8*int_sz + total_size_h1*dbl_sz) + 8*int_sz,'bof') == -1)
                    display('Error seeking Pressure file');
                    finish;
                end
                                    
                p = fread(fid,total_size_h1,'double');
                fclose(fid);
            
                p3_loc = zeros(l_ni_h1,l_nj_h1,l_nk_h1);
                %transform flat to 3D array
                for i = 1:l_ni_h1
                    for j = 1:l_nj_h1
                        for k = 1:l_nk_h1
                            p3_loc(i,j,k) = p(l_ni_h1*l_nj_h1*(k-1) + l_ni_h1*(j-1) + i);                         
                        end
                    end
                end
            
                %combining data from diff CPUs in X-dir
                if pi==1 
                    p_glob = p3_loc(2:end-1,2:end-1,2:end-1);                
                else
                    p_glob = [p_glob;p3_loc(2:end-1,2:end-1,2:end-1)];                
                end
            end
            %combining data from diff CPUs in Y-dir
            if pj==1 
                pp_glob = p_glob;
            else
                pp_glob = cat(2, pp_glob, p_glob);
            end
        end
        %combining data from diff CPUs in Z-dir
        if pk==1 
            P = pp_glob;
        else
            P = cat(3, P, pp_glob);
        end
    end
    clear p; clear p3_loc; clear p_glob; clear pp_glob;
    
    % display P
    figure;
    %h=pcolor(squeeze(Xp(:,:,z_slice)), squeeze(Yp(:,:,z_slice)), squeeze(P(:,:,z_slice)));       
    %set( h, 'EdgeColor', 'none');
    %colorbar;    
    [C,h] = contour(squeeze(Xp(:,:,z_slice)), squeeze(Yp(:,:,z_slice)), squeeze(P(:,:,z_slice)));       
    set(h,'ShowText','on','TextStep',get(h,'LevelStep'));   
    %colormap cool;        
    hold on;
    %bathymetry    
    line(Xp(:,1,z_slice),Yp(:,1,z_slice), 'Color', 'black', 'LineWidth', 2);
    axis('image');
    title('Fluid pressure $$P$$: longitudal $$z$$-slice','Interpreter','LaTex','FontSize',16);
end

% LOAD VELOCITY
if(plot_velocity_x_slice || plot_velocity_y_slice || plot_velocity_z_slice)
    
    [U, V, W] = load_binary_velocity(directory,timestep, delta_ts/save_timestep_period, npx,npy,npz, l_ni_h2,l_nj_h2,l_nk_h2);
     
    % PLOT VELOCITY: Z-SLICE
    if(plot_velocity_z_slice)            
        figure;    
        xz = squeeze(Xp(:,:,z_slice)); yz = squeeze(Yp(:,:,z_slice)); 
        uz = squeeze(U(:,:,z_slice));  vz = squeeze(V(:,:,z_slice)); wz = squeeze(W(:,:,z_slice));         
        uz = uz - mean(mean(uz));
        vz = vz - mean(mean(vz));
        wz = wz - mean(mean(wz));
        cav = abs(curl(yz,xz,vz,uz));
        f = cav;
        h=pcolor(xz, yz, f); set( h, 'EdgeColor', 'none');
        %contour(xz, yz, f, 50); 
        colorbar;
        %colormap bone;caxis([200 max(max(f))-min(min(f))])
        hold on;
        quiver(xz, yz, uz, vz, 'Color', [0 0 0]);         
        %streamslice(yz,xz, uz,vz, 10); 
        %bathymetry    
        line(xz(:,1),yz(:,1), 'Color', 'black', 'LineWidth', 2);
        axis('image');
        set(gca, 'FontSize', 14);   
        xlabel('$$X$$','Interpreter','LaTex','FontSize',16);     
        ylabel('$$Y$$','Interpreter','LaTex','FontSize',16);   
        title('Fluid velocity $$(u,v)$$: longitudal $$z$$-slice','Interpreter','LaTex','FontSize',16);
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