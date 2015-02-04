function show_function_modify(max_pull, r_param, order)
    [x y] = meshgrid(-10:0.1:10,-10:0.1:10);
    z = max_pull*exp(-r_param*(x.^2+y.^2).^(order/2));
    figure
    surfc(x,y,z);
    colorbar
    figure
    contour(x,y,z,15);
    colorbar
end