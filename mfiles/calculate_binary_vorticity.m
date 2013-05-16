function omega_1 = calculate_binary_vorticity(x,y,z,v,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates vorticity in the longitudinal (x) direction using binary
% output from cns-koltakov. x,y,z are the grid points (cell centers) 
% and v,w are the cartesian velocities at the cell centers.
%
% Bobby Arthur
% May 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate metrics
X_xi = 0.5*( x(3:end,2:end-1,2:end-1) - x(1:end-2,2:end-1,2:end-1) );
Y_xi = 0.5*( y(3:end,2:end-1,2:end-1) - y(1:end-2,2:end-1,2:end-1) );
Z_xi = 0.5*( z(3:end,2:end-1,2:end-1) - z(1:end-2,2:end-1,2:end-1) );

X_et = 0.5*( x(2:end-1,3:end,2:end-1) - x(2:end-1,1:end-2,2:end-1) );
Y_et = 0.5*( y(2:end-1,3:end,2:end-1) - y(2:end-1,1:end-2,2:end-1) );
Z_et = 0.5*( z(2:end-1,3:end,2:end-1) - z(2:end-1,1:end-2,2:end-1) );

X_zt = 0.5*( x(2:end-1,2:end-1,3:end) - x(2:end-1,2:end-1,1:end-2) );
Y_zt = 0.5*( y(2:end-1,2:end-1,3:end) - y(2:end-1,2:end-1,1:end-2) );
Z_zt = 0.5*( z(2:end-1,2:end-1,3:end) - z(2:end-1,2:end-1,1:end-2) );

J = ( X_xi.*Y_et.*Z_zt + X_et.*Y_zt.*Z_xi + X_zt.*Y_xi.*Z_et ...
    - X_zt.*Y_et.*Z_xi - X_et.*Y_xi.*Z_zt - X_xi.*Y_zt.*Z_et );

% XI_X = ( Y_et.*Z_zt - Y_zt.*Z_et );
% ET_X = ( Y_zt.*Z_xi - Y_xi.*Z_zt );
% ZT_X = ( Y_xi.*Z_et - Y_et.*Z_xi );
XI_Y = ( Z_et.*X_zt - Z_zt.*X_et );
ET_Y = ( Z_zt.*X_xi - Z_xi.*X_zt );
ZT_Y = ( Z_xi.*X_et - Z_et.*X_xi );
XI_Z = ( X_et.*Y_zt - X_zt.*Y_et );
ET_Z = ( X_zt.*Y_xi - X_xi.*Y_zt );
ZT_Z = ( X_xi.*Y_et - X_et.*Y_xi );


%Calculate vorticity
omega_1 = 1/J.* ...
    ( XI_Y.*( w(3:end,2:end-1,2:end-1) - w(1:end-2,2:end-1,2:end-1) )/2 ...
    + ET_Y.*( w(2:end-1,3:end,2:end-1) - w(2:end-1,1:end-2,2:end-1) )/2 ...
    + ZT_Y.*( w(2:end-1,2:end-1,3:end) - w(2:end-1,2:end-1,1:end-2) )/2 ...
    - XI_Z.*( v(3:end,2:end-1,2:end-1) - v(1:end-2,2:end-1,2:end-1) )/2 ...
    - ET_Z.*( v(2:end-1,3:end,2:end-1) - v(2:end-1,1:end-2,2:end-1) )/2 ...
    - ZT_Z.*( v(2:end-1,2:end-1,3:end) - v(2:end-1,2:end-1,1:end-2) )/2 );

end




