%==========================================================================
% Loads grid nodes from binary output
%==========================================================================
function [Nx, Ny, Nz, npx, npy, npz, Nt, period, ...
    dt, x_length, y_length, z_length] = load_binary_parameters(directory)
    
filename = [directory 'parameters.0'];   
    
str = sprintf('Reading Parameters...');  disp(str);
  
fid = fopen([filename],'r'); 

ints = fread(fid,8,'int');
Nx = ints(1);
Ny = ints(2);
Nz = ints(3);
npx = ints(4);
npy = ints(5);
npz = ints(6);
Nt = ints(7);
period = ints(8);

doubles = fread(fid,4,'double');
dt = doubles(1);
x_length = doubles(2);
y_length = doubles(3);
z_length = doubles(4);
fclose(fid);

end


                