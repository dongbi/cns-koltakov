%==========================================================================
% Loads Pressure from binary output
%
% if (delta_ts >0), an averaging is performed over delta_ts frames
%==========================================================================

function [P] = load_binary_pressure(directory, ts_start,delta_ts, npx,npy,npz, l_ni_h1,l_nj_h1,l_nk_h1)
    filename = [directory 'pressure.'];
    total_size_h1 = l_ni_h1 * l_nj_h1 * l_nk_h1;     % flat array_3d size w/ 2 halo cells
    
    int_sz = 4;
    dbl_sz = 8;
    
    for pk = 1:npz
        for pj = 1:npy
            for pi = 1:npx 
                proc = (pi-1)*npz*npy + (pj-1)*npz + pk-1;
                str = sprintf('Proc = %d',proc);  disp(str); 
                fid = fopen([filename num2str(proc)],'r');
                
                p3_loc = zeros(l_ni_h1,l_nj_h1,l_nk_h1);                     
                
                % averaging if (delta_ts > 1)
                for frame = ts_start:ts_start+delta_ts                
                    if(fseek(fid,frame*(8*int_sz + total_size_h1*dbl_sz) + 8*int_sz,'bof') == -1)
                        display('Error seeking binary file: Pressure');
                        quit;
                    end               
                                    
                    p_flat = fread(fid,total_size_h1,'double');                               

                    %transform flat to 3D array
                    for i = 1:l_ni_h1
                        for j = 1:l_nj_h1
                            for k = 1:l_nk_h1
                                p3_loc(i,j,k) = p3_loc(i,j,k) + p_flat(l_ni_h1*l_nj_h1*(k-1) + l_ni_h1*(j-1) + i); 
                            end
                        end
                    end                                        
                end %for:frame
                
                % divide by the number of frames if averaging is required
                if(delta_ts>0)
                    p3_loc(:,:,:) = p3_loc(:,:,:) ./ delta_ts;
                end
                
                fclose(fid);
                
                %combining data from diff CPUs in X-dir
                if pi==1 
                    p_glob = p3_loc(2:end-1,2:end-1,2:end-1);
                else
                   p_glob  = [p_glob;p3_loc(2:end-1,2:end-1,2:end-1)];
                end
            end
            %combining data from diff CPUs in Y-dir
            if pj==1 
                pp_glob = p_glob;
            else
                pp_glob = cat(2, pp_glob, p_glob);
            end
        end
        %combining data from diff CPUs in Z-dirfseek(fid,8*int_sz,'bof');                
        if pk==1 
            P = pp_glob;
        else
            P = cat(3, P, pp_glob);
        end
    end
    clear p_flat; clear p3_loc; clear p_glob; clear pp_glob;