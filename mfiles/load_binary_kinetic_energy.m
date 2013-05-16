%==========================================================================
% Loads kinetic_energy from binary output
%==========================================================================

function E_k = load_binary_kinetic_energy(directory, ts_start)
    filename = [directory 'kinetic_energy.0'];
    dbl_sz = 8;
    
    fid = fopen(filename,'r');                   
%      keyboard;         
    if(fseek(fid,ts_start*dbl_sz,'bof') == -1)
        display('Error seeking binary file: Kinetic Energy');
        return;
    end               

    E_k = fread(fid,1,'double'); 

    fclose(fid);
end
                
                