%==========================================================================
% Loads Velocities from binary output
%
% if (delta_ts >0), an averaging is performed over delta_ts frames
%==========================================================================

function [U, V, W] = load_binary_velocity(directory, ts_start,delta_ts, npx,npy,npz, l_ni_h2,l_nj_h2,l_nk_h2)
    filename = [directory 'velocity.'];
    total_size_h2 = l_ni_h2 * l_nj_h2 * l_nk_h2;     % flat array_3d size w/ 2 halo cells
    
    int_sz = 4;
    dbl_sz = 8;
    
    str = sprintf('Reading Velocity...');  disp(str);
    for pk = 1:npz
        for pj = 1:npy
            for pi = 1:npx 
                proc = (pi-1)*npz*npy + (pj-1)*npz + pk-1;
%                 str = sprintf('Proc = %d',proc);  disp(str); 
                fid = fopen([filename num2str(proc)],'r');
                
                u = zeros(l_ni_h2,l_nj_h2,l_nk_h2); v = u; w = u;                     
                
                % averaging if (delta_ts > 1)
                for frame = ts_start:ts_start+delta_ts                
%                     if(fseek(fid,frame*(8*int_sz + 3*total_size_h2*dbl_sz) + 8*int_sz,'bof') == -1)
                      if(fseek(fid,frame*(3*total_size_h2*dbl_sz),'bof') == -1)
                        display('Error seeking binary file: Velocity');
                        return;
                    end               
                                    
                    uvw = fread(fid,3*total_size_h2,'double');                               

                    %transform flat to 3D array
                    for i = 1:l_ni_h2
                        for j = 1:l_nj_h2
                            for k = 1:l_nk_h2
                                u(i,j,k) = u(i,j,k) + uvw(                l_ni_h2*l_nj_h2*(k-1) + l_ni_h2*(j-1) + i); 
                                v(i,j,k) = v(i,j,k) + uvw(total_size_h2  +l_ni_h2*l_nj_h2*(k-1) + l_ni_h2*(j-1) + i); 
                                w(i,j,k) = w(i,j,k) + uvw(2*total_size_h2+l_ni_h2*l_nj_h2*(k-1) + l_ni_h2*(j-1) + i); 
                            end
                        end
                    end                                        
                end %for:frame
                
                % divide by the number of frames if averaging is required
                if(delta_ts>0)
                    u(:,:,:) = u(:,:,:) ./ delta_ts;
                    v(:,:,:) = v(:,:,:) ./ delta_ts;
                    w(:,:,:) = w(:,:,:) ./ delta_ts;
                end
                
                fclose(fid);
                
                %combining data from diff CPUs in X-dir
                if pi==1 
                    u_glob = u(3:end-2,3:end-2,3:end-2);
                    v_glob = v(3:end-2,3:end-2,3:end-2);
                    w_glob = w(3:end-2,3:end-2,3:end-2);
                else
                   u_glob  = [u_glob;u(3:end-2,3:end-2,3:end-2)];
                   v_glob  = [v_glob;v(3:end-2,3:end-2,3:end-2)];
                   w_glob  = [w_glob;w(3:end-2,3:end-2,3:end-2)];
                end
            end
            %combining data from diff CPUs in Y-dir
            if pj==1 
                uu_glob = u_glob; vv_glob = v_glob; ww_glob = w_glob;
            else
                uu_glob = cat(2, uu_glob, u_glob); vv_glob = cat(2, vv_glob, v_glob); ww_glob = cat(2, ww_glob, w_glob);
            end
        end
        %combining data from diff CPUs in Z-dirfseek(fid,8*int_sz,'bof');                
        if pk==1 
            U = uu_glob; V = vv_glob; W = ww_glob;
        else
            U = cat(3, U, uu_glob); V = cat(3, V, vv_glob); W = cat(3, W, ww_glob);
        end
    end
    clear uvw;
    clear u; clear v; clear w;
    clear u_glob; clear v_glob; clear w_glob;
    clear uu_glob; clear vv_glob; clear ww_glob;               