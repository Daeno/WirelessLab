[x y] = meshgrid(-10:0.1:10,-10:0.1:10);
z = 0.4*exp(-0.3*(x.^2+y.^2));
figure
surfc(x,y,z);
colorbar
figure
contour(x,y,z,15);
colorbar