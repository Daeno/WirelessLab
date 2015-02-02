function [x,y] = init_pos(n,rangex,rangey)
    xrange=rangex(2)-rangex(1);
    yrange=rangey(2)-rangey(1);
    x=rand(1,n)*xrange+rangex(1);
    y=rand(1,n)*yrange+rangey(1);
end