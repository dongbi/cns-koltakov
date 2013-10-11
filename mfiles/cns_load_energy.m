% Calculates and saves energy timeseries for binary CNS output
clear all; clc; close all;

directory = '/output/';

%Load options
filename = 'energy.mat';
timestep_initial = 25;
timestep_final = 10000;
delta_ts = 0; %always 0

%Parameters
[Nx, Ny, Nz, npx, npy, npz, Nt, save_timestep_period, ...
    dt, x_length, y_length, z_length] = load_binary_parameters(directory);
g = 9.81;
nu = 10^-6;
kappa = 10^-6;

% local grid dims
% halo <- 2
l_ni_h2 = Nx/npx + 4;
l_nj_h2 = Ny/npy + 4;
l_nk_h2 = Nz/npz + 4;

%Adjust timestep
t = (timestep_initial:save_timestep_period:timestep_final);
Nsteps = length(t); 
timesteps = (t - save_timestep_period)/save_timestep_period;
t = t*dt;

%Grid
[xh,yh,zh] = load_binary_grid_w_halo(directory,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2,1); %with halo=1
[XI_X,ET_X,ZT_X,XI_Y,ET_Y,ZT_Y,XI_Z,ET_Z,ZT_Z,J] = calculate_binary_metrics(xh,yh,zh);

%Load energy variables
Eb = zeros(1,Nsteps); Ep = Eb; Ek = Eb; phi_d = Eb; phi_z = Eb; phi_i = Eb; epsilon_tot = Eb;
Ek1 = Eb; Ek2 = Eb; Ek3 = Eb;
F_Eb = Eb; F_Ep = Ep; F_Ek = Eb; %not used 

for n=1:Nsteps
    str = sprintf(['n = ',num2str(n),', t = ',num2str(t(n))]); disp(str);
    tn = timesteps(n);
    
    %Ep,Eb,Ek
    [Eb(n),Ep(n),phi_d(n),F_Eb(n),F_Ep(n)] = load_binary_potential_energy(directory,tn);
    [Ek(n),F_Ek(n),epsilon_tot(n)] = load_binary_kinetic_energy(directory,tn);

    %phi_i
    [rho_n] = load_binary_density(directory,tn,...
      delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
    phi_i(n) = kappa*g* ( sum(sum(squeeze(rho_n(:,:,end)).*squeeze(ZT_Z(:,:,end)))) ...
                        -sum(sum(squeeze(rho_n(:,:,  1)).*squeeze(ZT_Z(:,:,  1)))) );

    %phi_z
    [~,~,w_n] = load_binary_velocity(directory,tn,...
      delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
    phi_z(n) = g*sum(sum(sum(rho_n.*w_n.*J)));

    clear 'rho_n' 'w_n';
end

%Adjust values
Eb = Eb-Eb(1);
Ep = Ep-Ep(1);
Ea = Ep - Eb;
Et = Ep + Ek;

%Take time derivatives
dEtdt = (Et(3:end)-Et(1:end-2))/(2*dt*save_timestep_period);
dEpdt = (Ep(3:end)-Ep(1:end-2))/(2*dt*save_timestep_period);
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*dt*save_timestep_period);
dEadt = (Ea(3:end)-Ea(1:end-2))/(2*dt*save_timestep_period);
dEkdt = (Ek(3:end)-Ek(1:end-2))/(2*dt*save_timestep_period);
t_dt = t(2:end-1);

%Save
save(filename,'Et','Ep','Eb','Ea','Ek','dEtdt','dEpdt','dEbdt','dEadt', ...
    'dEkdt','phi_d','epsilon_tot','phi_z','phi_i','t','t_dt');
display('Complete');
