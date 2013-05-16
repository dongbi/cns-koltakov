function lambda_2 = calculate_binary_lambda_2(x,y,z,u,v,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates lambda_2 using binary output from cns-koltakov. x,y,z are the 
% grid points (cell centers) and u,v,w are the cartesian velocities at the 
% cell centers. lambda_2 is the median eigenvalue of the tensor 
% (S_ik*S_kj + Omega_ik*Omega_kj) where S and Omega are the symmetric and
% antisymmetric parts of the strain rate tensor du_i/dx_j.
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

XI_X = ( Y_et.*Z_zt - Y_zt.*Z_et );
ET_X = ( Y_zt.*Z_xi - Y_xi.*Z_zt );
ZT_X = ( Y_xi.*Z_et - Y_et.*Z_xi );
XI_Y = ( Z_et.*X_zt - Z_zt.*X_et );
ET_Y = ( Z_zt.*X_xi - Z_xi.*X_zt );
ZT_Y = ( Z_xi.*X_et - Z_et.*X_xi );
XI_Z = ( X_et.*Y_zt - X_zt.*Y_et );
ET_Z = ( X_zt.*Y_xi - X_xi.*Y_zt );
ZT_Z = ( X_xi.*Y_et - X_et.*Y_xi );

%Calculate components of strain rate s = du_i/dx_j
s_11 = 1/J.* ...
     ( XI_X.*( u(3:end,2:end-1,2:end-1) - u(1:end-2,2:end-1,2:end-1) )/2 ...
     + ET_X.*( u(2:end-1,3:end,2:end-1) - u(2:end-1,1:end-2,2:end-1) )/2 ...
     + ZT_X.*( u(2:end-1,2:end-1,3:end) - u(2:end-1,2:end-1,1:end-2) )/2 );
 
s_12 = 1/J.* ...
     ( XI_Y.*( u(3:end,2:end-1,2:end-1) - u(1:end-2,2:end-1,2:end-1) )/2 ...
     + ET_Y.*( u(2:end-1,3:end,2:end-1) - u(2:end-1,1:end-2,2:end-1) )/2 ...
     + ZT_Y.*( u(2:end-1,2:end-1,3:end) - u(2:end-1,2:end-1,1:end-2) )/2 );
 
s_13 = 1/J.* ...
     ( XI_Z.*( u(3:end,2:end-1,2:end-1) - u(1:end-2,2:end-1,2:end-1) )/2 ...
     + ET_Z.*( u(2:end-1,3:end,2:end-1) - u(2:end-1,1:end-2,2:end-1) )/2 ...
     + ZT_Z.*( u(2:end-1,2:end-1,3:end) - u(2:end-1,2:end-1,1:end-2) )/2 );
 
s_21 = 1/J.* ...
     ( XI_X.*( v(3:end,2:end-1,2:end-1) - v(1:end-2,2:end-1,2:end-1) )/2 ...
     + ET_X.*( v(2:end-1,3:end,2:end-1) - v(2:end-1,1:end-2,2:end-1) )/2 ...
     + ZT_X.*( v(2:end-1,2:end-1,3:end) - v(2:end-1,2:end-1,1:end-2) )/2 );
 
s_22 = 1/J.* ...
     ( XI_Y.*( v(3:end,2:end-1,2:end-1) - v(1:end-2,2:end-1,2:end-1) )/2 ...
     + ET_Y.*( v(2:end-1,3:end,2:end-1) - v(2:end-1,1:end-2,2:end-1) )/2 ...
     + ZT_Y.*( v(2:end-1,2:end-1,3:end) - v(2:end-1,2:end-1,1:end-2) )/2 );
 
s_23 = 1/J.* ...
     ( XI_Z.*( v(3:end,2:end-1,2:end-1) - v(1:end-2,2:end-1,2:end-1) )/2 ...
     + ET_Z.*( v(2:end-1,3:end,2:end-1) - v(2:end-1,1:end-2,2:end-1) )/2 ...
     + ZT_Z.*( v(2:end-1,2:end-1,3:end) - v(2:end-1,2:end-1,1:end-2) )/2 );

s_31 = 1/J.* ...
     ( XI_X.*( w(3:end,2:end-1,2:end-1) - w(1:end-2,2:end-1,2:end-1) )/2 ...
     + ET_X.*( w(2:end-1,3:end,2:end-1) - w(2:end-1,1:end-2,2:end-1) )/2 ...
     + ZT_X.*( w(2:end-1,2:end-1,3:end) - w(2:end-1,2:end-1,1:end-2) )/2 );
 
s_32 = 1/J.* ...
     ( XI_Y.*( w(3:end,2:end-1,2:end-1) - w(1:end-2,2:end-1,2:end-1) )/2 ...
     + ET_Y.*( w(2:end-1,3:end,2:end-1) - w(2:end-1,1:end-2,2:end-1) )/2 ...
     + ZT_Y.*( w(2:end-1,2:end-1,3:end) - w(2:end-1,2:end-1,1:end-2) )/2 );
 
s_33 = 1/J.* ...
     ( XI_Z.*( w(3:end,2:end-1,2:end-1) - w(1:end-2,2:end-1,2:end-1) )/2 ...
     + ET_Z.*( w(2:end-1,3:end,2:end-1) - w(2:end-1,1:end-2,2:end-1) )/2 ...
     + ZT_Z.*( w(2:end-1,2:end-1,3:end) - w(2:end-1,2:end-1,1:end-2) )/2 );

%Calculate lambda_2
Nx = size(x,1); Ny = size(x,2); Nz = size(x,3);
lambda_2 = zeros(Nx-2,Ny-2,Nz-2);

for i=1:Nx-2
    for j=1:Ny-2
        for k=1:Nz-2
            %Build S and Omega
            s = [s_11(i,j,k), s_12(i,j,k), s_13(i,j,k); ...
                 s_21(i,j,k), s_22(i,j,k), s_23(i,j,k); ...
                 s_31(i,j,k), s_32(i,j,k), s_33(i,j,k)];
            S = 0.5*(s+s');
            Omega = 0.5*(s-s');
            
            %Calculate lambda_2
            lambda = eig(S*S + Omega*Omega);
            lambda_2(i,j,k) = lambda(2);
        end
    end
end
end




