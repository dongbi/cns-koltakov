% Saves frames for an isosurface movie of CNS simulation data loaded from binary 
% output files
clear all; clc; close all;

directory = '/usr/var/tmp/barthur/1109466/output/';

% PLOTTING OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timestep_initial = 50;
timestep_final = 14800;
delta_ts = 0;

rho_iso = 0;
omega_1_iso = 1.25;

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

%Adjust timestep
t = (timestep_initial:save_timestep_period:timestep_final);
timestep = (t - save_timestep_period)/save_timestep_period;

%grid
[x,y,z] = load_binary_grid(directory,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);

% ISOSURFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Bottom
z_slice = round(Nz/2);
xz = squeeze(x(:,:,z_slice)); yz = squeeze(y(:,:,z_slice)); 
depth = abs(squeeze(z(:,:,1)));

%Figure properties
h1=figure(1);
set(h1,'visible','off');
set(h1,'Position',[100 100 1500 1000]);
set(h1,'PaperPositionMode','auto');
axis off;
hold on;

for n=timestep
    [rho] = load_binary_density(directory,n,...
           delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
    cla;
    surf(xz,yz,-depth-.005,'FaceColor','k','EdgeColor','none'); 
    pat = patch(isosurface(x,y,z,rho,rho_iso,rho));
    set(pat,'FaceColor','interp','EdgeColor','none');
    set(pat,'FaceColor','r');
    
%     [u,v,w] = load_binary_velocity(directory,n,...
%        delta_ts/save_timestep_period, npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
%     omega_1 = calculate_binary_vorticity(x,y,z,v,w);
%     pat_pos = patch(isosurface( ...
%                  x(2:end-1,2:end-1,2:end-1), ...
%                  y(2:end-1,2:end-1,2:end-1), ...
%                  z(2:end-1,2:end-1,2:end-1), ...
%                  omega_1,omega_1_iso,omega_1) );
%     set(pat_pos,'FaceColor','b','EdgeColor','none');
%     pat_neg = patch(isosurface( ...
%                  x(2:end-1,2:end-1,2:end-1), ...
%                  y(2:end-1,2:end-1,2:end-1), ...
%                  z(2:end-1,2:end-1,2:end-1), ...
%                  omega_1,-omega_1_iso,omega_1) );
%     set(pat_neg,'FaceColor','g','EdgeColor','none');
    
    daspect([1,1,1]);
    view([0.04, -0.04, 0.04]); axis tight;
    camlight ;
    lighting gouraud;
    
    saveas(isoplot,['isoplot_',num2str(n+1),'fig']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('Complete');
