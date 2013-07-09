%==========================================================================
% Loads kinetic_energy from binary output
%==========================================================================

function [Ek,F_Ek,epsilon] = load_binary_kinetic_energy(directory, ts_start)
    filename = [directory 'kinetic_energy.0'];
    dbl_sz = 8;
    
    fid = fopen(filename,'r');                         
    if(fseek(fid,ts_start*(3*dbl_sz),'bof') == -1)
        display('Error seeking binary file: Kinetic Energy');
        return;
    end               

    kinetic_energy = fread(fid,3,'double'); 
    Ek = kinetic_energy(1);
    F_Ek = kinetic_energy(2);
    epsilon = kinetic_energy(3);
    fclose(fid);
end
                
                