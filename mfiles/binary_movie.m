load '/home/barthur/zang/christine_12deg_2/rho.mat'

%%
figure;
hold on;
axis equal;
for n=1:254
    cla;
    pcolor(x(:,y_slice,:),z(:,y_slice,:),rho(:,:,n));
    shading flat;
    drawnow;
    pause(0.1);
end




        
    