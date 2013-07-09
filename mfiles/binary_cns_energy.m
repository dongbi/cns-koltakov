% Calculates and saves energy timeseries for binary CNS output
clear all; clc; close all;

directory = '/home/barthur/zang/3D_test/';
% directory = '/usr/var/tmp/barthur/1118351/output/';
timestep_initial = 1;
timestep_final = 300;
load_timestep_period = 1;
delta_ts = 0;

%Parameters
[Nx, Ny, Nz, npx, npy, npz, Nt, save_timestep_period, ...
    dt, x_length, y_length, z_length] = load_binary_parameters(directory);
g = 9.81;
nu = 10^-6;

% local grid dims
% halo <- 2
l_ni_h2 = Nx/npx + 4;
l_nj_h2 = Ny/npy + 4;
l_nk_h2 = Nz/npz + 4;

%Adjust timestep
t = (timestep_initial:save_timestep_period:timestep_final);
Nsteps = length(t); 
t_load = (load_timestep_period:load_timestep_period:timestep_final);
Nload = length(t_load);
timesteps = (t - save_timestep_period)/save_timestep_period;
t = t*dt;
t_load = t_load*dt;

%Grid
[x,y,z] = load_binary_grid(directory,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2); %without halo
yl = y(1,:,:); zl = z(1,:,:); %left boundary grid points
[xh,yh,zh] = load_binary_grid_w_halo(directory,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2,1); %with halo=1
[XI_X,ET_X,ZT_X,XI_Y,ET_Y,ZT_Y,XI_Z,ET_Z,ZT_Z,J] = calculate_binary_metrics(xh,yh,zh);
Al = (1/J(1,:,:)) / (x_length/Nx); %y-z area of left boundary grid cells

%Load energy variables
Eb = zeros(1,Nsteps); Ep = Eb; Ek = Eb; 
F_Eb = zeros(1,Nsteps); F_Ep = F_Eb; phi_d = F_Eb; F_Ek = F_Eb; epsilon = F_Eb;
epsilon_tot = zeros(1,Nload); epsilon_tke = epsilon_tot;

for n=1:Nsteps
    str = sprintf(['n = ',num2str(n),', t = ',num2str(t(n))]); disp(str);
    tn = timesteps(n);
    
    %Ep,Eb,Ek
    [Eb(n),Ep(n),phi_d(n),F_Eb(n),F_Ep(n)] = load_binary_potential_energy(directory,tn);
    [Ek(n),F_Ek(n),epsilon(n)] = load_binary_kinetic_energy(directory,tn);
    
    if(mod(n,load_timestep_period)==0)
%         %Boundary PE flux
%         ul = 0.03*cos(pi/z_length*zl)*sin(2*pi/10*t(n));
%         [rho] = load_binary_density(directory,tn,...
%             delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
%         rhol = rho(1,:,:);
%         clear 'rho';
%     
%         F_pe(n/load_timestep_period) = g*sum(sum(rhol.*zl.*ul.*Al));
%     
%         %Boundary KE flux
%         F_ke(n/load_timestep_period) = sum(sum(rhol.*(ul.^3).*Al));
    
        %Dissipation
        [uh,vh,wh] = load_binary_velocity_w_halo(directory, tn,delta_ts/save_timestep_period, npx,npy,npz, l_ni_h2,l_nj_h2,l_nk_h2, 1);
        epsilon_tot(n/load_timestep_period) = calculate_binary_dissipation(uh,vh,wh,nu,XI_X,ET_X,ZT_X,XI_Y,ET_Y,ZT_Y,XI_Z,ET_Z,ZT_Z,J);
%         epsilon_tke(n/load_timestep_period) = calculate_binary_dissipation(uh-repmat(mean(uh(:,2:end-1,:),2),[1 Ny+2 1]),...
%                                                   vh-repmat(mean(vh(:,2:end-1,:),2),[1 Ny+2 1]),...
%                                                   wh-repmat(mean(wh(:,2:end-1,:),2),[1 Ny+2 1]),...
%                                                   nu,XI_X,ET_X,ZT_X,XI_Y,ET_Y,ZT_Y,XI_Z,ET_Z,ZT_Z,J);
        clear 'uh' 'vh' 'wh';
    end
end

%Adjust values
Eb = Eb-Eb(1);
Ep = Ep-Ep(1);
Ea = Ep - Eb;
Et = Ep + Ek;

%Take time derivatives
dEtdt = (Et(3:end)-Et(1:end-2))/(2*dt);
dEpdt = (Ep(3:end)-Ep(1:end-2))/(2*dt);
dEbdt = (Eb(3:end)-Eb(1:end-2))/(2*dt);
dEadt = (Ea(3:end)-Ea(1:end-2))/(2*dt);
dEkdt = (Ek(3:end)-Ek(1:end-2))/(2*dt);
t_dt = t(2:end-1);

%Save
save('energy1.mat','Et','Ep','Eb','Ea','Ek','dEtdt','dEpdt','dEbdt','dEadt', ...
     'dEkdt','F_Eb','F_Ep','phi_d','F_Ek','epsilon_tot','epsilon_tke','t','t_dt','t_load');
display('Complete');
