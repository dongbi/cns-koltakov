function omega_1 = calculate_binary_vorticity(XI_Y,ET_Y,ZT_Y,XI_Z,ET_Z,ZT_Z,J,v,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates vorticity in the longitudinal (x) direction using binary
% output from cns-koltakov. The metric quantities are calculated using 
% the function calculate_binary_metrics.m and v,w are the cartesian 
% velocities at the cell centers including 1 halo cell.
%
% Bobby Arthur
% May 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate vorticity
omega_1 = 1/J.* ...
    ( XI_Y.*( w(3:end,2:end-1,2:end-1) - w(1:end-2,2:end-1,2:end-1) )/2 ...
    + ET_Y.*( w(2:end-1,3:end,2:end-1) - w(2:end-1,1:end-2,2:end-1) )/2 ...
    + ZT_Y.*( w(2:end-1,2:end-1,3:end) - w(2:end-1,2:end-1,1:end-2) )/2 ...
    - XI_Z.*( v(3:end,2:end-1,2:end-1) - v(1:end-2,2:end-1,2:end-1) )/2 ...
    - ET_Z.*( v(2:end-1,3:end,2:end-1) - v(2:end-1,1:end-2,2:end-1) )/2 ...
    - ZT_Z.*( v(2:end-1,2:end-1,3:end) - v(2:end-1,2:end-1,1:end-2) )/2 );

end




