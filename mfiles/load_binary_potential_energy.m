%==========================================================================
% Loads potential_energy from binary output
%==========================================================================

function [E_b,E_p] = load_binary_potential_energy(directory, ts_start)
    filename = [directory 'potential_energy.0'];
    dbl_sz = 8;
    
    fid = fopen(filename,'r');                   
%      keyboard;         
    if(fseek(fid,ts_start*(2*dbl_sz),'bof') == -1)
        display('Error seeking binary file: Potential Energy');
        return;
    end               

    potential_energy = fread(fid,2,'double'); 
    E_b = potential_energy(1);
    E_p = potential_energy(2);

    fclose(fid);
end
                
                