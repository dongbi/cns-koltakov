% Displays CNS simulation data from binary output files
clear all; clc; close all;
% directory = '/home/barthur/zang/3D_test/';
% directory = '/home/barthur/zang/2D_test/';
directory = '/home/barthur/zang/christine_12deg_2/';

% PLOTTING OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timestep = 12700; %timestep to plot
delta_ts = 0; %averaging

display_grid = 0;
display_density = 0;
display_velocity = 0;
display_scalar = 0;
display_pressure = 0;
display_density_isosurface = 1;
    rho_iso = -.005;
    
x_loc = 1.85; %if 0, west boundary
y_loc = 0; %if 0, centerline

plot_xz = 1; %x-z plot
plot_yz = 1; %y-z plot

plot_pcolor = 1;
plot_quiver = 0;
show_grid_lines = 0; %false = shading flat

% OTHER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Nx, Ny, Nz, npx, npy, npz, Nt, save_timestep_period, ...
    dt, x_length, y_length, z_length] = load_binary_parameters(directory);

% local grid dims
% halo <- 1
l_ni_h1 = Nx/npx + 2;
l_nj_h1 = Ny/npy + 2;
l_nk_h1 = Nz/npz + 2;
% halo <- 2
l_ni_h2 = Nx/npx + 4;
l_nj_h2 = Ny/npy + 4;
l_nk_h2 = Nz/npz + 4;

total_size_h1 = l_ni_h1 * l_nj_h1 * l_nk_h1;  % flat array_3d size w/ 1 halo cell
total_size_h2 = l_ni_h2 * l_nj_h2 * l_nk_h2;  % flat array_3d size w/ 2 halo cells
int_sz = 4;
dbl_sz = 8;

%adjust timestep
timestep = timestep - save_timestep_period;
timestep = timestep/save_timestep_period;

% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,y,z] = load_binary_grid(directory,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);

if(x_loc)
    x_slice = find(squeeze(x(:,1,1))>=x_loc,1,'first');
else
    x_slice = 1;
end

if(y_loc)
    y_slice = find(squeeze(y(1,:,1))>=y_loc,1,'first');
else
    y_slice = round(Ny/2); 
end

if(plot_xz)
   x_xz = squeeze(x(:,y_slice,:)); 
   z_xz = squeeze(z(:,y_slice,:)); 
end

if(plot_yz)
    y_yz = squeeze(y(x_slice,:,:));
    z_yz = squeeze(z(x_slice,:,:));
end

if(display_grid)  
    if(plot_xz)
        figure;
        plot(x_xz(:,:), z_xz(:,:), 'k.');     
        axis equal;
        axis([0 x_length -z_length 0]);
        xlabel('x [m]');
        ylabel('z [m]');
        title('Computational grid, x-z slice');
    end
    
    if(plot_yz)
        figure;
        plot(y_yz(:,:), z_yz(:,:), 'k.');         
        axis equal;
        axis([0 y_length -z_length 0]);
        xlabel('y [m]');
        ylabel('z [m]');
        title('Computational grid, y-z slice');
    end
end

% DENSITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_density)
   [rho] = load_binary_density(directory,timestep,...
       delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
   
   if(plot_xz)
       rho_xz = squeeze(rho(:,y_slice,:));
       figure;
       pcolor(x_xz,z_xz,rho_xz);
       colorbar;
       axis equal;
       axis([0 x_length -z_length 0]);
       xlabel('x [m]');
       ylabel('z [m]');
       title('Density \Delta\rho/\rho_0, x-z slice');
       if(~show_grid_lines)
           shading flat;
       end
       if (plot_yz)
           hold on;
           plot([x_loc x_loc],[z(x_slice,1,1) z(x_slice,1,end)],'k-','LineWidth',2);
           hold off;
       end
   end
   
   if(plot_yz)
       rho_yz = squeeze(rho(x_slice,:,:));
       figure;
       pcolor(y_yz,z_yz,rho_yz);
       colorbar;
       axis equal;
       axis([0 y_length -z_length 0]);
       xlabel('y [m]');
       ylabel('z [m]');
       title('Density \Delta\rho/\rho_0, y-z slice');
       if(~show_grid_lines)
           shading flat;
       end
   end
end

% VELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_velocity)
    [u,v,w] = load_binary_velocity(directory,timestep,...
        delta_ts/save_timestep_period, npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
    
    if(plot_xz)
        u_xz = squeeze(u(:,y_slice,:));
        v_xz = squeeze(v(:,y_slice,:));
        w_xz = squeeze(w(:,y_slice,:));
        figure;
        pcolor(x_xz,z_xz,u_xz);
        colorbar;
%         quiver(x_xz,z_xz,u_xz,w_xz);
        axis equal;
        axis([0 x_length -z_length 0]);
        xlabel('x [m]');
        ylabel('z [m]');
        title('Velocity u, x-z slice');
        if(~show_grid_lines)
            shading flat;
        end
    end
    
    if(plot_yz)
        u_yz = squeeze(u(x_slice,:,:));
        v_yz = squeeze(v(x_slice,:,:));
        w_yz = squeeze(w(x_slice,:,:));
        figure;
        pcolor(y_yz,z_yz,u_yz);
        colorbar;
        axis equal;
        axis([0 y_length -z_length 0]);
        xlabel('y [m]');
        ylabel('z [m]');
        title('Velocity u, y-z slice');
        if(~show_grid_lines)
            shading flat;
        end
    end  
end

% SCALAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_scalar)
   [phi] = load_binary_scalar(directory,timestep,...
       delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);

   if(plot_xz)
       phi_xz = squeeze(phi(:,y_slice,:));
       figure;
       pcolor(x_xz,z_xz,phi_xz);
       colorbar;
       caxis([0 1]);
       axis equal;
       axis([0 x_length -z_length 0]);
       xlabel('x [m]');
       ylabel('z [m]');
       title('Scalar \phi, x-z slice');
       if(~show_grid_lines)
           shading flat;
       end
   end
   
   if(plot_yz)
       phi_yz = squeeze(phi(x_slice,:,:));
       figure;
       pcolor(y_yz,z_yz,phi_yz);
       colorbar;
       caxis([0 1]);
       axis equal;
       axis([0 y_length -z_length 0]);
       xlabel('y [m]');
       ylabel('z [m]');
       title('Scalar \phi, y-z slice');
       if(~show_grid_lines)
           shading flat;
       end
   end
end

% PRESSURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_pressure)    
    [P] = load_binary_pressure(directory, timestep,delta_ts, npx,npy,npz, l_ni_h1,l_nj_h1,l_nk_h1);

    if(plot_xz)
       P_xz = squeeze(P(:,y_slice,:));
       figure;
       pcolor(x_xz,z_xz,P_xz);
       colorbar;
       axis equal;
       axis([0 x_length -z_length 0]);
       xlabel('x [m]');
       ylabel('z [m]');
       title('Pressure, x-z slice');
       if(~show_grid_lines)
           shading flat;
       end
   end
   
   if(plot_yz)
       phi_yz = squeeze(P(x_slice,:,:));
       figure;
       pcolor(y_yz,z_yz,P_yz);
       colorbar;
       axis equal;
       axis([0 y_length -z_length 0]);
       xlabel('y [m]');
       ylabel('z [m]');
       title('Pressure, y-z slice');
       if(~show_grid_lines)
           shading flat;
       end
   end
end

% DENSITY ISOSURFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_density_isosurface)
    z_slice = round(Nz/2);
    xz = squeeze(x(:,:,z_slice)); yz = squeeze(y(:,:,z_slice)); 
    depth = abs(squeeze(z(:,:,1)));

    [rho] = load_binary_density(directory,timestep,...
           delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);

    figure;
    surf(xz,yz,-depth-.005,'FaceColor','k','EdgeColor','none'); 
    hold on;
    pat = patch(isosurface(x,y,z,rho,0,rho));
    set(pat,'FaceColor','interp','EdgeColor','none');
    daspect([1,1,1])
    view([0.04, -0.04, 0.04]); axis tight
    camlight 
    lighting gouraud
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    title(['Isosurface of \Delta\rho/\rho_0=',num2str(rho_iso)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf, 'PaperPositionMode', 'auto');  