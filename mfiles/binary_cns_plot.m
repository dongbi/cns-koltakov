% Plots CNS simulation data from .mat file created with
% binary_cns_load.m
clear all; clc; 

% load /home/barthur/Desktop/1109466_Eb.mat
% load /home/barthur/Desktop/freeslip.mat
% load /home/barthur/Desktop/gdu11550.mat
load /home/barthur/zang/GRC/midslice-no-3.mat
% load ./test.mat
% load /home/barthur/Documents/MATLAB/research/zang_binary_out/test2.mat

% PLOTTING OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display_grid = 0;
display_density = 1;
display_velocity = 0;
display_scalar = 0;
display_pressure = 0;
display_isosurface = 0;
display_density_isosurface = 0;
    rho_iso = 0;
display_omega_1_isosurface = 0;
    omega_1_iso = 1.25;
display_lambda_2_isosurface = 0;
    lambda_2_iso = -0.1;
display_Q_criterion_isosurface = 0;
    Q_iso = 1;
display_potential_energy = 0;
    
x_loc = 0; %if 0, west boundary
y_loc = 0; %if 0, centerline

slices = 1;
plot_xz = 1; %x-z plot
plot_yz = 0; %y-z plot

plot_quiver = 0;
    iskip = 3;
    kskip = 3;
show_grid_lines = 0; %false = shading flat

% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(slices)
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
       u_xz = squeeze(u(:,y_slice,:));
       v_xz = squeeze(v(:,y_slice,:));
       w_xz = squeeze(w(:,y_slice,:));
    end

    if(plot_yz)
        y_yz = squeeze(y(x_slice,:,:));
        z_yz = squeeze(z(x_slice,:,:));
        u_yz = squeeze(u(x_slice,:,:));
        v_yz = squeeze(v(x_slice,:,:));
        w_yz = squeeze(w(x_slice,:,:));
    end
end

if(display_grid)  
    if(plot_xz)
        grid_fig_xz = figure;
        hold all;
        plot(x_xz(:,:), z_xz(:,:), 'k.');     
        axis equal;
        axis([0 x_length -z_length 0]);
        xlabel('x [m]');
        ylabel('z [m]');
        title('Computational grid, x-z slice');
    end
    
    if(plot_yz)
        grid_fig_yz = figure;
        hold all;
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
   
   if(plot_xz)
       rho_xz = squeeze(rho(:,y_slice,:));
       rho_fig_xz = figure;
       hold all;
       pcolor(x_xz,z_xz,rho_xz);
       colorbar;
       if(plot_quiver)
            quiver(x_xz(1:iskip:end,1:kskip:end),z_xz(1:iskip:end,1:kskip:end),u_xz(1:iskip:end,1:kskip:end),w_xz(1:iskip:end,1:kskip:end),'k');
       end
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
       rho_fig_yz = figure;
       hold all;
       pcolor(y_yz,z_yz,rho_yz);
       colorbar;
       if(plot_quiver)
            quiver(y_yz(1:iskip:end,1:kskip:end),z_yz(1:iskip:end,1:kskip:end),v_yz(1:iskip:end,1:kskip:end),w_yz(1:iskip:end,1:kskip:end),'k');
       end
       axis image;
%        axis([0 y_length -z_length 0]);
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
    
    if(plot_xz)
        u_xz = squeeze(u(:,y_slice,:));
        v_xz = squeeze(v(:,y_slice,:));
        w_xz = squeeze(w(:,y_slice,:));
        velocity_fig_xz = figure;
        hold all;
        pcolor(x_xz,z_xz,u_xz);
        colorbar;
        if(plot_quiver)
            quiver(x_xz(1:iskip:end,1:kskip:end),z_xz(1:iskip:end,1:kskip:end),u_xz(1:iskip:end,1:kskip:end),w_xz(1:iskip:end,1:kskip:end),'k');
        end
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
        velocity_fig_yz = figure;
        hold all;
        pcolor(y_yz,z_yz,u_yz);
        colorbar;
        if(plot_quiver)
            quiver(y_yz(1:iskip:end,1:kskip:end),z_yz(1:iskip:end,1:kskip:end),v_yz(1:iskip:end,1:kskip:end),w_yz(1:iskip:end,1:kskip:end),'k');
        end
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

   if(plot_xz)
       phi_xz = squeeze(phi(:,y_slice,:));
       phi_fig_xz = figure;
       hold all;
       pcolor(x_xz,z_xz,phi_xz);
       colorbar;
       caxis([0 1]);
       if(plot_quiver)
           quiver(x_xz(1:iskip:end,1:kskip:end),z_xz(1:iskip:end,1:kskip:end),u_xz(1:iskip:end,1:kskip:end),w_xz(1:iskip:end,1:kskip:end),'k');
       end
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
       phi_fig_yz = figure;
       hold all;
       pcolor(y_yz,z_yz,phi_yz);
       colorbar;
       caxis([0 1]);
       if(plot_quiver)
           quiver(y_yz(1:iskip:end,1:kskip:end),z_yz(1:iskip:end,1:kskip:end),v_yz(1:iskip:end,1:kskip:end),w_yz(1:iskip:end,1:kskip:end),'k');
       end
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

    if(plot_xz)
       p_xz = squeeze(p(:,y_slice,:));
       p_fig_xz = figure;
       hold all;
       pcolor(x_xz,z_xz,p_xz);
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
       phi_yz = squeeze(p(x_slice,:,:));
       p_fig_yz = figure;
       hold all;
       pcolor(y_yz,z_yz,p_yz);
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

% ISOSURFACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_isosurface)
    z_slice = round(Nz/2);
    xz = squeeze(x(:,:,z_slice)); yz = squeeze(y(:,:,z_slice)); 
    depth = abs(squeeze(z(:,:,1)));

    isoplot = figure;
    hold all;
    bot = surf(xz,yz,-depth,'FaceColor','k','EdgeColor','none'); 
    axis off;
    
    if(display_density_isosurface)
        pat = patch(isosurface(x,y,z,rho,rho_iso,rho));
        set(pat,'FaceColor','interp','EdgeColor','none');
        set(pat,'FaceColor','r');
    end
    
    if(display_omega_1_isosurface)
        omega_1 = calculate_binary_vorticity(x,y,z,v,w);
        pat_pos = patch(isosurface( ...
                     x(2:end-1,2:end-1,2:end-1), ...
                     y(2:end-1,2:end-1,2:end-1), ...
                     z(2:end-1,2:end-1,2:end-1), ...
                     omega_1,omega_1_iso,omega_1) );
        set(pat_pos,'FaceColor','b','EdgeColor','none');
        pat_neg = patch(isosurface( ...
                     x(2:end-1,2:end-1,2:end-1), ...
                     y(2:end-1,2:end-1,2:end-1), ...
                     z(2:end-1,2:end-1,2:end-1), ...
                     omega_1,-omega_1_iso,omega_1) );
        set(pat_neg,'FaceColor','g','EdgeColor','none');
    end
    
    if(display_lambda_2_isosurface)
        lambda_2 = calculate_binary_lambda_2(x,y,z,u,v,w);
        pat_lambda = patch(isosurface( ...
             x(2:end-1,2:end-1,2:end-1), ...
             y(2:end-1,2:end-1,2:end-1), ...
             z(2:end-1,2:end-1,2:end-1), ...
             lambda_2,lambda_2_iso,lambda_2) );
         set(pat_lambda,'FaceColor','b','EdgeColor','none');
    end
    
    if(display_Q_criterion_isosurface)
        Q = calculate_binary_Q_criterion(x,y,z,u,v,w);
        pat_Q = patch(isosurface( ...
             x(2:end-1,2:end-1,2:end-1), ...
             y(2:end-1,2:end-1,2:end-1), ...
             z(2:end-1,2:end-1,2:end-1), ...
             Q,Q_iso,Q) );
         set(pat_Q,'FaceColor','b','EdgeColor','none');
    end
    
    daspect([1,1,1])
    view([0.04, -0.04, 0.04]); axis tight
    camlight 
    lighting gouraud
%     xlabel('x [m]');
%     ylabel('y [m]');
%     zlabel('z [m]');
%     title(['Isosurface of \Delta\rho/\rho_0=',num2str(rho_iso)]);
end

% POTENTIAL ENERGY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_potential_energy)
    PE_fig = figure(2);
    hold all;
    Eb = Eb - Eb(1);
    Ep = Ep - Ep(1);
    Ea = Ep - Eb;
    Et = Ep + Ek;
    
%     Ebmin = -findpeaks(-Eb(1:1500));
%     Ebmin = Ebmin(1);
%     Epmin = -findpeaks(-Ep(1:1500));
%     Ebamp = (Eb(1) - Ebmin)/2;
%     Epamp = (Ep(1) - Epmin)/2;
%     Eb = Eb - Ebamp*cos(2*pi/10*t) + Ebamp;
%     Ep = Ep - Epamp*cos(2*pi/10*t) + Epamp;
    
%     t = t/dt;
%     plot(t(1:14000),Eb(1:14000),'b',t(1:14000),Ep(1:14000),'k',t(1:14000),Ea(1:14000),'r');
%     hold on;
%     plot([3700 3700],[-1 1],'k:');
%     plot([5500 5500],[-1 1],'k:');
%     plot([7600 7600],[-1 1],'k:');
%     plot([9600 9600],[-1 1],'k:');
%     plot([11700 11700],[-1 1],'k:');
%     plot([13800 13800],[-1 1],'k:');
    plot(t,Ep,'k--',t,Eb,'b',t,Ea,'r',t,Ek,'k:',t,Et,'k');
    xlabel('t [s]');
    ylabel('E-E_0');
    title('Time-series of Potential and Kinetic Energy');
    legend('E_p','E_b','E_a','E_k','E_T');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Complete');





        
    