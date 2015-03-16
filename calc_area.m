function [ctr,aMap] = calc_area(x,y,z,ranges,steps,points)
    l_y = size(z,1);
    l_x = size(z,2);
    aMap = zeros(l_y, l_x);
    for i = 1:l_y
        for j= 1:l_x
            if(z(i,j) < 0.1)
                aMap(i,j) = 1;
            end
        end
    end
    xmin = ranges(1,1);
    xmax = ranges(1,2);
    ymin = ranges(2,1);
    ymax = ranges(2,2);
    d_x = [-1,-1,-1,0,0,0,1,1,1];
    d_y = [-1,0,1,-1,0,1,-1,0,1];
    for k = 1:size(points,1)
        x_index = ceil((points(k,1) - xmin)/steps(1));
        y_index = ceil((points(k,2) - ymin)/steps(2));
        for i=1:9
            n_x = x_index + d_x(i);
            n_y = y_index + d_y(i);
            if(n_x > 0 && n_y > 0 && n_x <= l_x && n_y <= l_y)
                aMap(n_y,n_x) = 1;
            end
        end
    end
    ctr = 0;
    for i = 1:l_y
        for j = 1:l_x
            if(aMap(i,j) == 1)
                ctr = ctr + 1;
            end
        end
    end
    
end