%==========================================================================
% Loads potential_energy from binary output
%==========================================================================

function [Eb,Ep,phi_d,F_Eb,F_Ep] = load_binary_potential_energy(directory, ts_start)
    filename = [directory 'potential_energy.0'];
    dbl_sz = 8;
    
    fid = fopen(filename,'r');                          
    if(fseek(fid,ts_start*(5*dbl_sz),'bof') == -1)
        display('Error seeking binary file: Potential Energy');
        return;
    end               

    potential_energy = fread(fid,5,'double'); 
    Eb = potential_energy(1);
    Ep = potential_energy(2);
    phi_d = potential_energy(3);
    F_Eb = potential_energy(4);
    F_Ep = potential_energy(5);

    fclose(fid);
end
                
                