function [N, count_n,x,y,z] = firefly_simple(instr,x,y,z,p_x,p_y)
%% parameter initialization
if nargin<1,   instr=[30 200 1];     end
n=instr(1);  MaxGeneration=instr(2);  draw = instr(3);

range=[-10 10 ; -10 10]; % range=[xmin xmax ymin ymax];
N = 0;
threshold = -1;

% constants for reduction of intensity
[ma, na] = meshgrid(-0.6:0.2:0.6, -0.6:0.2:0.6);
ra = ma.^2 + na.^2;
a = 5.*power(4, -ra);

%% algorithm-related parameters
alpha           =   0.2 * ones(1,instr(1));     % Randomness, bigger the more random
step            =   0.2;                        % Maximum step a node can move in one iteration
rpl_range       =   0.8;                        % Repulsive force range
gamma           =   1;                          % Distance decay constant
delta           =   0.7;                        % 
% function modifing parameters
global max_pull r_param order      
max_pull        =   0.3;
r_param         =   0.3;
order           =   3;
show_function_modify(max_pull, r_param, order);
%% target function generation
if(nargin < 4)
[x,y,z] = randfunc_g(range,100);
end

figure(1);
surfc(x,y,z);
%% deciding initial point
if(nargin < 6)
    [xn, yn]=init_pos(n,range(1,:),range(2,:));
else
    xn = p_x;
    yn = p_y;
end

holdTime = zeros(size(xn, 2));
maxHoldTime = 3;
xd = [];
yd = [];

lowArea = (z < -1);
sumLowArea = sum(sum(lowArea))

figure(2);
count_n = 0;

%% Main Loop
zn = interp2(x, y, z, xn, yn);
while(count_n < MaxGeneration)
    count_n         =   count_n + 1;
    xo              =   xn;   
    yo              =   yn;
    zo              =   zn;%interp2(x, y, z, xo, yo);
    max_zo          =   max(zo);

    % Compute the undiscovered area; 
    ztemp = z;
    ni = size(xn,2) + size(xd,2);
    gridrange = size(x);
    xtemp = [xn xd];
    ytemp = [yn yd];
    [I, J] = findindex(xtemp, ytemp, range, step);
    for i=1:ni,
        for m = I(i)-3:I(i)+3,
            for n = J(i)-3:J(i)+3,
                if m > 1 && m <= gridrange(1, 1) && n > 1 && n <= gridrange(1, 2),
                    ztemp(n,m) = ztemp(n,m) + a(m-I(i)+4, n-J(i)+4);
                end
            end
        end
    end
    lowArea = (ztemp < -1);
    sumLowArea = sum(sum(lowArea));
    N(count_n)      =   sumLowArea;

%% Move all fireflies to the better locations
    [xn,yn,xd,yd,holdTime]  =   ffa_move(xo, yo, zo, step, alpha, gamma, range, threshold, rpl_range, xd, yd, holdTime, maxHoldTime);
    zn                      =   interp2(x, y, z, xn, yn);
    alpha                   =   newalpha(alpha,delta,threshold,zn);
%% Drawing utility
    if(draw)
        axis equal;
        contour(x,y,ztemp,15);
        hold on;
        plot(xo(zo>threshold),yo(zo>threshold),'.','markersize',10,'markerfacecolor','r');
        plot(xo(zo<=threshold),yo(zo<=threshold),'g.','markersize',10,'markerfacecolor','r');
        quiver(xo,yo,xn-xo,yn-yo,'AutoScale','off');
        drawnow;
        hold off;
    end
    
    if(max_zo < threshold) break;  end
end % end of while
xd
yd
end
%% ----- All subfunctions are listed here ---------

%% Move all fireflies toward brighter ones
function [xn,yn,xd,yd,holdTime] = ffa_move(xo, yo, zo, step, alpha, gamma, range, threshold, rpl_range, xd, yd, holdTime, maxHoldTime)
    ni=size(yo,2); nj=size(yo,2);
    xn = xo;
    yn = yo;
    %% Brightness modification
    %  based on r^2
    zo = modifing_target_f(xo, yo, zo, xd, yd);
    min_zo = min(zo);
    for i=1:ni,
        %% whether to stop when fulfilled the criteria
        if zo(i) < threshold, 
            holdTime(i) = holdTime(i) + 1;
        %        continue; 
        elseif holdTime(i) >= maxHoldTime,
            holdTime(i) = 0;
            xd = [xd xn(i)];
            yd = [yd yn(i)];
        else
            holdTime(i) = 0;
        end
        nd = size(xd,2);
        for d=1:nd,
            r   =   sqrt((xo(i)-xd(d))^2+(yo(i)-yd(d))^2);
            if r < rpl_range,
                xn(i)= xn(i) + (xo(i) - xd(d));
                yn(i)= yn(i) + (yo(i) - yd(d));
            end
        end
        for j=1:nj,
            r   =   sqrt((xo(i)-xo(j))^2+(yo(i)-yo(j))^2);
            %% simple repulsive force node i if j is in the range of 
            %if r < rpl_range,
            %    xn(i)= xn(i) + (xo(i) - xo(j));
            %    yn(i)= yn(i) + (yo(i) - yo(j));
            %    continue;
            %end
            %% Attractivness based on Ligntness(or z) (z larger, lightness deemer)
            if  zo(i) > zo(j) && zo(j) < threshold,     % Brighter and smaller than 
                beta0   =   min_zo/zo(j);               % Normalize Brightness
                beta    =   beta0 * exp(-gamma*r.^2);   % exponential decay with r^2
                xn(i)   =   xn(i) +(xo(j) - xo(i))*beta;
                yn(i)   =   yn(i) +(yo(j) - yo(i))*beta;
            end
        end % end for j
        r_g = rand * step;
        theta = rand * 2 * pi;
        xn(i) = xn(i) + alpha(i)*r_g*cos(theta);
        yn(i) = yn(i) + alpha(i)*r_g*sin(theta);
    end % end for i

    dx = xn - xo;
    dy = yn - yo;
    r = dx.*dx + dy.*dy;

    for idx = 1:numel(dx)
       if(r(idx) > step^2)
            dx(idx) = dx(idx) / r(idx).^0.5 * step;
            dy(idx) = dy(idx) / r(idx).^0.5 * step;
       end
    end

    xn = dx + xo;
    yn = dy + yo;
    [xn,yn] = findrange(xn,yn,range);                   % confine xn,yn inside the range
end
%% Reduce the randomness during iterations
% maximum alpha = 3
function alpha=newalpha(alpha,delta,threshold,zn)
    for i=1:numel(zn)
        if(zn(i) <= threshold)
            alpha(i)=alpha(i)*delta;
        else
            alpha(i) = alpha(i) / delta;
        end
        if(alpha(i) > 3)
            alpha(i) = 3;
        end
    end
end
%% function modification
function zn = modifing_target_f(xo, yo, zo, xd, yd)
    zn = zo;
    ni = size(zo,2);
    nj = ni;
    nd = size(xd,2);
    global max_pull r_param order
    for i = 1:ni,
        for j = 1:nj,
            r   =  sqrt((xo(i)-xo(j))^2+(yo(i)-yo(j))^2);
            zn(i) = zn(i) + max_pull*exp(-r_param*r^order)*zo(j);
        end
        for j = i:nd,
            r   =  sqrt((xo(i)-xd(j))^2+(yo(i)-yd(j))^2);
            zn(i) = zn(i) + max_pull*exp(-r_param*r^order)*9;
        end
    end
end

function [I, J] = findindex(p_x, p_y, range, step)
    I = ceil((p_x-range(1, 1))/step);
    J = ceil((p_y-range(2, 1))/step);
end