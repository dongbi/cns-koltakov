function Q = calculate_binary_Q_criterion(XI_X,ET_X,ZT_X,XI_Y,ET_Y,ZT_Y,...
                                          XI_Z,ET_Z,ZT_Z,J,u,v,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates Q-criterion using binary output from cns-koltakov. The metric 
% quantities are calculated using the function calculate_binary_metrics.m 
% and u,v,w are the cartesian velocities at the cell centers including 1 
% halo cell. 
% 
% Q=1/2*(||Omega||^2 - ||S||^2) 
%
% where ||A|| = sqrt(tr(AA')) and S and Omega are the symmetric and
% antisymmetric parts of the strain rate tensor du_i/dx_j.
%
% Bobby Arthur
% May 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%Calculate Q-criterion
Nx = size(u,1); Ny = size(u,2); Nz = size(u,3);
Q = zeros(Nx-2,Ny-2,Nz-2);

for i=1:Nx-2
    for j=1:Ny-2
        for k=1:Nz-2
            %Build S and Omega
            s = [s_11(i,j,k), s_12(i,j,k), s_13(i,j,k); ...
                 s_21(i,j,k), s_22(i,j,k), s_23(i,j,k); ...
                 s_31(i,j,k), s_32(i,j,k), s_33(i,j,k)];
            S = 0.5*(s+s');
            Omega = 0.5*(s-s');
            
            %Calculate Q
            Q(i,j,k) = 0.5*( abs(trace(Omega*Omega')) - abs(trace(S*S')) );
        end
    end
end

end


