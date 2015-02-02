function coverage = cover_percent(x, y, z, p_x, p_y, range, threshold ,node_range)
%% A function that caculate the percentage of low point that has been found
%  caculate the coverage in following criteria
%  
zz = z < threshold;
low_count = sum(zz);
n = size(p_x,2);
table = zeros(size(x));
%% naive method
for i = 1:n
    dx = x - p_x(i);
    dy = y - p_y(i);
    r = dx.^2 + dy.^2;
    table = table | ((r<node_range)&&(zz));
end
coverage = sum(table)/low_count;
end