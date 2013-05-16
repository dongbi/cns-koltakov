% Loads CNS simulation data from binary output files to .mat file
clear all; clc; close all;

% directory = '/usr/var/tmp/barthur/1101642/output/';
% directory = '/usr/var/tmp/barthur/1102291/output/';
% directory = '/usr/var/tmp/barthur/1109177/output/';
% directory = '/usr/var/tmp/barthur/1109466/output/';
% directory = '/home/barthur/Desktop/';
directory = '/home/barthur/zang/3D_test/';

% LOAD OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'test2.mat';
timestep_initial = 2;
timestep_final = 300; 
delta_ts = 0; %averaging

load_grid = 0;
load_density = 0;
load_velocity = 0;
load_scalar = 0;
load_pressure = 0;
load_potential_energy = 1;

lateral_average = 1;

% LOAD STATIC PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose variable to save
save_vars = '''Nx'',''Ny'',''Nz'',''Nt'',''save_timestep_period'',''dt''';
save_vars = [save_vars,',''t'',''x_length'',''y_length'',''z_length'''];
if(load_grid)
    save_vars = [save_vars,',''x'',''y'',''z'''];
end
if(load_density)
    save_vars = [save_vars,',''rho'''];
end
if(load_velocity)
    save_vars = [save_vars,',''u'',''v'',''w'''];
end
if(load_scalar)
    save_vars = [save_vars,',''phi'''];
end
if(load_pressure)
    save_vars = [save_vars,',''p'''];
end
if(load_potential_energy)
    save_vars = [save_vars,',''Eb'',''Ep'',''Ek'''];
end
if(lateral_average)
    save_vars = [save_vars,',''lateral_average'''];
end

%Parameters
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
Nsteps = length(t);
timesteps = (t - save_timestep_period)/save_timestep_period;
t = t*dt;

%Grid
if(load_grid)
    [x,y,z] = load_binary_grid(directory,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
    if(lateral_average)
        x = x(:,1,:);
        y = y(:,1,:);
        z = z(:,1,:);
    end
end

% LOAD TIME VARIABLE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize storage variables
if(lateral_average)
    rho = zeros(Nx,Nz,Nsteps);
    u = rho; v = rho; w = rho; phi = rho; p = rho;
else
    rho = zeros(Nx,Ny,Nz,Nsteps);
    u = rho; v = rho; w = rho; phi = rho; p = rho;
end
Eb = zeros(1,Nsteps); Ep = Eb; Ek = Eb;

for n=1:Nsteps
    str = sprintf(['n = ',num2str(n),', t = ',num2str(t(n))]); disp(str);
    tn = timesteps(n);
    
    %Density
    if(load_density)
        [rho_n] = load_binary_density(directory,tn,...
           delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
        if(lateral_average)
            rho(:,:,n) = mean(rho_n,2);
        else
            rho(:,:,:,n) = rho_n;
        end
        clear 'rho_n';
    end
    
    %Velocity
    if(load_velocity)
        [u_n,v_n,w_n] = load_binary_velocity(directory,tn,...
            delta_ts/save_timestep_period, npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
        if(lateral_average)
            u(:,:,n) = mean(u_n,2);
            v(:,:,n) = mean(v_n,2);
            w(:,:,n) = mean(w_n,2);
        else
            u(:,:,:,n) = u_n;
            v(:,:,:,n) = v_n;
            w(:,:,:,n) = w_n;
        end
        clear 'u_n' 'v_n' 'w_n';
    end
    
    %Scalar
    if(load_scalar)
        [phi_n] = load_binary_scalar(directory,tn,...
            delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
        if(lateral_average)
            phi(:,:,n) = mean(phi_n,2);
        else
            phi(:,:,:,n) = phi_n;
        end
        clear 'phi_n';
    end
    
    %Pressure
    if(load_pressure)    
        [p_n] = load_binary_pressure(directory,tn,...
            delta_ts, npx,npy,npz, l_ni_h1,l_nj_h1,l_nk_h1);
        if(lateral_average)
            p(:,:,n) = mean(p_n,2);
        else
            p(:,:,:,n) = p_n;
        end
        clear 'p_n';
    end
    
    %Potential Energy
    if(load_potential_energy)
        [Eb(n),Ep(n)] = load_binary_potential_energy(directory,tn);
        Ek(n) = load_binary_kinetic_energy(directory,tn);
    end
end

%Squeeze out extra dimensions
rho = squeeze(rho);
u = squeeze(u);
v = squeeze(v);
w = squeeze(w);
phi = squeeze(phi);
p = squeeze(p);

%Save 
eval(['save(''',filename,''',',save_vars,')']);
display('Complete');

