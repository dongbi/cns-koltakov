%==========================================================================
% Loads Velocities from binary output and calculates RMS velocity
%
% if (delta_ts >0), an averaging is performed over delta_ts frames
%==========================================================================

function [Up, Vp, Wp, UpVp, UpWp, VpWp, TKE] = load_binary_rms_velocity(directory, ts_start,delta_ts, U_ave, V_ave, W_ave, npx,npy,npz, l_ni_h2,l_nj_h2,l_nk_h2)
    filename = [directory 'velocity.'];
    total_size_h2 = l_ni_h2 * l_nj_h2 * l_nk_h2;     % flat array_3d size w/ 2 halo cells
    
    l_ni = l_ni_h2 - 4;
    l_nj = l_nj_h2 - 4;
    l_nk = l_nk_h2 - 4;
    
    int_sz = 4;
    dbl_sz = 8;
    
    for pk = 1:npz
        for pj = 1:npy
            for pi = 1:npx 
                proc = (pi-1)*npz*npy + (pj-1)*npz + pk-1;
                str = sprintf('Proc = %d',proc);  disp(str); 
                fid = fopen([filename num2str(proc)],'r');
                
                u = zeros(l_ni_h2,l_nj_h2,l_nk_h2); v=u; w=u; uv=u; uw=u; vw=u; tke=u;                     
                
                % averaging if (delta_ts > 1)
                for frame = ts_start:ts_start+delta_ts                
                    if(fseek(fid,frame*(8*int_sz + 3*total_size_h2*dbl_sz) + 8*int_sz,'bof') == -1)
                        display('Error seeking binary file: Velocity');
                        quit;
                    end               
                                    
                    uvw = fread(fid,3*total_size_h2,'double');                               

                    %transform flat to 3D array (operate on internal cells)
                    for i = 3:l_ni_h2-2
                        for j = 3:l_nj_h2-2
                            for k = 3:l_nk_h2-2
				u_diff = uvw(                l_ni_h2*l_nj_h2*(k-1) + l_ni_h2*(j-1) + i) - U_ave((pi-1)*l_ni+i-2,(pj-1)*l_nj+j-2,(pk-1)*l_nk+k-2);
                                v_diff = uvw(total_size_h2  +l_ni_h2*l_nj_h2*(k-1) + l_ni_h2*(j-1) + i) - V_ave((pi-1)*l_ni+i-2,(pj-1)*l_nj+j-2,(pk-1)*l_nk+k-2);
				w_diff = uvw(2*total_size_h2+l_ni_h2*l_nj_h2*(k-1) + l_ni_h2*(j-1) + i) - W_ave((pi-1)*l_ni+i-2,(pj-1)*l_nj+j-2,(pk-1)*l_nk+k-2);
				u(i,j,k) = u(i,j,k) + u_diff^2; 
                                v(i,j,k) = v(i,j,k) + v_diff^2; 
                                w(i,j,k) = w(i,j,k) + w_diff^2;
				uv(i,j,k) = uv(i,j,k) + u_diff*v_diff;
				uw(i,j,k) = uw(i,j,k) + u_diff*w_diff;
				vw(i,j,k) = vw(i,j,k) + v_diff*w_diff;
				tke(i,j,k) = tke(i,j,k) + u_diff^2 + v_diff^2 + w_diff^2; 
                            end
                        end
                    end                                        
                end %for:frame
                
                % divide by the number of frames if averaging is required
                if(delta_ts>0)
                    u(:,:,:) = sqrt(u(:,:,:) ./ delta_ts);
                    v(:,:,:) = sqrt(v(:,:,:) ./ delta_ts);
                    w(:,:,:) = sqrt(w(:,:,:) ./ delta_ts);
		    uv(:,:,:) = -uv(:,:,:) ./ delta_ts;
                    uw(:,:,:) = -uw(:,:,:) ./ delta_ts;
                    vw(:,:,:) = -vw(:,:,:) ./ delta_ts;
		    tke(:,:,:) = tke(:,:,:) ./ (2*delta_ts);
                end
                
                fclose(fid);
                
                %combining data from diff CPUs in X-dir
                if pi==1 
                    u_glob = u(3:end-2,3:end-2,3:end-2);
                    v_glob = v(3:end-2,3:end-2,3:end-2);
                    w_glob = w(3:end-2,3:end-2,3:end-2);
		    uv_glob = uv(3:end-2,3:end-2,3:end-2);
                    uw_glob = uw(3:end-2,3:end-2,3:end-2);
                    vw_glob = vw(3:end-2,3:end-2,3:end-2);
                    tke_glob = tke(3:end-2,3:end-2,3:end-2);
                else
                   u_glob  = [u_glob;u(3:end-2,3:end-2,3:end-2)];
                   v_glob  = [v_glob;v(3:end-2,3:end-2,3:end-2)];
                   w_glob  = [w_glob;w(3:end-2,3:end-2,3:end-2)];
		   uv_glob  = [uv_glob;uv(3:end-2,3:end-2,3:end-2)];
                   uw_glob  = [uw_glob;uw(3:end-2,3:end-2,3:end-2)];
                   vw_glob  = [vw_glob;vw(3:end-2,3:end-2,3:end-2)];
                   tke_glob  = [tke_glob;tke(3:end-2,3:end-2,3:end-2)];
                end
            end
            %combining data from diff CPUs in Y-dir
            if pj==1 
                uu_glob = u_glob; vv_glob = v_glob; ww_glob = w_glob; 
		uvuv_glob = uv_glob; uwuw_glob = uw_glob; vwvw_glob = vw_glob; tketke_glob = tke_glob;
            else
                uu_glob = cat(2, uu_glob, u_glob); vv_glob = cat(2, vv_glob, v_glob); ww_glob = cat(2, ww_glob, w_glob);
            	uvuv_glob = cat(2, uvuv_glob, uv_glob); uwuw_glob = cat(2, uwuw_glob, uw_glob); vwvw_glob = cat(2, vwvw_glob, vw_glob);
		tketke_glob = cat(2, tketke_glob, tke_glob); 
	    end
        end
        %combining data from diff CPUs in Z-dirfseek(fid,8*int_sz,'bof');                
        if pk==1 
            Up = uu_glob; Vp = vv_glob; Wp = ww_glob; 
	    UpVp = uvuv_glob; UpWp = uwuw_glob; VpWp = vwvw_glob; 
            TKE = tketke_glob;
        else
            Up = cat(3, Up, uu_glob); Vp = cat(3, Vp, vv_glob); Wp = cat(3, Wp, ww_glob);
	    UpVp = cat(3, UpVp, uvuv_glob); UpWp = cat(3, UpWp, uwuw_glob); VpWp = cat(3, VpWp, vwvw_glob); 
            TKE = cat(3, TKE, tketke_glob); 
        end
    end
    clear uvw;
    clear u; clear v; clear w; clear uv; clear uw; clear vw; clear tke;
    clear u_glob; clear v_glob; clear w_glob; clear uv_glob; clear uw_glob; clear vw_glob; clear tke_glob;
    clear uu_glob; clear vv_glob; clear ww_glob; clear uvuv_glob; clear uwuw_glob; clear vwvw_glob; clear tketke_glob;
