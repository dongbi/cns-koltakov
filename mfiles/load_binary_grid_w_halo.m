%==========================================================================
% Loads grid nodes from binary output
%==========================================================================
function [Xp, Yp, Zp] = load_binary_grid_w_halo(directory, npx,npy,npz, l_ni_h2,l_nj_h2,l_nk_h2, hs)
    filename = [directory 'grid.'];
    total_size_h2 = l_ni_h2 * l_nj_h2 * l_nk_h2;     % flat array_3d size w/ 2 halo cells    
    int_sz = 4;   
    h = 2; %actual halo, hs is halo to save
    
    str = sprintf('Reading Grid...');  disp(str);
    for pk = 1:npz
        for pj = 1:npy
            for pi = 1:npx
                proc = (pi-1)*npz*npy + (pj-1)*npz + pk-1;
                fid = fopen([filename num2str(proc)],'r');   
%                 fseek(fid,8*int_sz,'bof');                
                g = fread(fid,3*total_size_h2,'double');
                fclose(fid);
                
                %transform flat to 3D array
                gx_loc = zeros(l_ni_h2,l_nj_h2,l_nk_h2); gy_loc = gx_loc; gz_loc = gy_loc;
                for i = 1:l_ni_h2
                    for j = 1:l_nj_h2
                        for k = 1:l_nk_h2
                            gx_loc(i, j, k) = g(                l_ni_h2*l_nj_h2*(k-1) + l_ni_h2*(j-1) + i); 
                            gy_loc(i, j, k) = g(total_size_h2  +l_ni_h2*l_nj_h2*(k-1) + l_ni_h2*(j-1) + i); 
                            gz_loc(i, j, k) = g(2*total_size_h2+l_ni_h2*l_nj_h2*(k-1) + l_ni_h2*(j-1) + i); 
                        end
                    end
                end
                
                %determine what part of array to keep
                first_x = 1+h;
                last_x = l_ni_h2-h;
                if pi==1
                    first_x = 1+hs;
                end
                if pi==npx
                    last_x = l_ni_h2-hs;
                end
                
                first_y = 1+h;
                last_y = l_nj_h2-h;
                if pj==1
                    first_y = 1+hs;
                end
                if pj==npy
                    last_y = l_nj_h2-hs;
                end
                
                first_z = 1+h;
                last_z = l_nk_h2-h; 
                if pk==1
                    first_z = 1+hs;
                end
                if pk==npz
                    last_z = l_nk_h2-hs;
                end 
                
                gx_loc = gx_loc(first_x:last_x,first_y:last_y,first_z:last_z);
                gy_loc = gy_loc(first_x:last_x,first_y:last_y,first_z:last_z);
                gz_loc = gz_loc(first_x:last_x,first_y:last_y,first_z:last_z);
                    
                %combining data from diff CPUs in X-dir
                if pi==1 
                    Gx = gx_loc;
                    Gy = gy_loc;
                    Gz = gz_loc;
                else
                    Gx = [Gx;gx_loc];
                    Gy = [Gy;gy_loc];
                    Gz = [Gz;gz_loc];
                end
            end
            %combining data from diff CPUs in Y-dir
            if pj==1 
                GGx = Gx; GGy = Gy; GGz = Gz;
            else
                GGx = cat(2, GGx, Gx); GGy = cat(2, GGy, Gy); GGz = cat(2, GGz, Gz);
            end
        end
        %combining data from diff CPUs in Z-dirfseek(fid,8*int_sz,'bof');                
        if pk==1 
            Xp = GGx; Yp = GGy; Zp = GGz;
        else
            Xp = cat(3, Xp, GGx); Yp = cat(3, Yp, GGy); Zp = cat(3, Zp, GGz);
        end
    end
    clear g;
    clear gx_loc; clear Gx; clear GGx;
    clear gy_loc; clear Gy; clear GGy;
    clear gz_loc; clear Gz; clear GGz;              