% Loads and saves CNS simulation data from binary output files
clear all; clc; close all;
% directory = '/home/barthur/zang/3D_test/';
% directory = '/home/barthur/zang/2D_test/';
directory = '/home/barthur/zang/christine_12deg_2/';
filename = '/home/barthur/zang/christine_12deg_2/rho.mat';
save_file = 1;

delta_ts = 0;

[Nx, Ny, Nz, npx, npy, npz, Nt, save_timestep_period, ...
    dt, x_length, y_length, z_length] = load_binary_parameters(directory);
Nt = 90000;

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

y_slice = round(Ny/2); 

%load grid
[x,y,z] = load_binary_grid(directory,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);

%load time dependent variables
rho = zeros(Nx,Nz,Nt/save_timestep_period);
u = rho; v = rho; w = rho; phi = rho; P = rho;

tstart = tic;
for n=0:Nt/save_timestep_period-1
    [rhon] = load_binary_density(directory,n,...
       delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
    rho(:,:,n+1) = rhon(:,y_slice,:);
%     [un,vn,wn] = load_binary_velocity(directory,n,...
%         delta_ts/save_timestep_period, npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
%     u(:,:,n+1) = un(:,y_slice,:);
%     w(:,:,n+1) = wn(:,y_slice,:);
%     [phi(:,:,n+1)] = load_binary_scalar(directory,n,...
%        delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
%     [P(:,:,n+1)] = load_binary_pressure(directory,n,delta_ts, npx,npy,npz, l_ni_h1,l_nj_h1,l_nk_h1);
    telapsed = toc(tstart);
    display(n);
    display(telapsed);
end

if(save_file)
   save(filename,'x','z','rho','y_slice');
end
    