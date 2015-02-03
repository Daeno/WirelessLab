function [N, count_n,x,y,z] = firefly_simple(instr,x,y,z,p_x,p_y)
%% parameter initialization
if nargin<1,   instr=[100 10000];     end
n=instr(1);  MaxGeneration=instr(2);  draw = instr(3);

range=[-10 10 ; -10 10]; % range=[xmin xmax ymin ymax];
N = 0;
threshold = -1;

%% algorithm-related parameters
alpha           =   0.2 * ones(1,instr(1));     % Randomness, bigger the more random
step            =   0.2;                        % Maximum step a node can move in one iteration
rpl_range       =   0.3;                        % Repulsive force range
gamma           =   1;                          % Distance decay constant
delta           =   0.7;                        % 
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

    N(count_n)      =   sum(zo < threshold);

%% Move all fireflies to the better locations
    [xn,yn]         =   ffa_move(xo, yo, zo, step, alpha, gamma, range, threshold, rpl_range);
    zn              =   interp2(x, y, z, xn, yn);
    alpha           =   newalpha(alpha,delta,threshold,zn);
%% Drawing utility
    if(draw)
        axis equal;
        contour(x,y,z,15);
        hold on;
        plot(xo,yo,'.','markersize',10,'markerfacecolor','g');
        quiver(xo,yo,xn-xo,yn-yo);
        drawnow;
        hold off;
    end
    
    if(max_zo < threshold) break;  end
end % end of while

end
%% ----- All subfunctions are listed here ---------

%% Move all fireflies toward brighter ones
function [xn,yn] = ffa_move(xo, yo, zo, step, alpha, gamma, range, threshold, rpl_range)
ni=size(yo,2); nj=size(yo,2);
xn = xo;
yn = yo;
zn = zo;
min_zo = min(zo);
for i=1:ni,
    %% whether to stop when fulfilled the criteria
    %if zo(i) < threshold, 
            %continue; 
    %end
    for j=1:nj,
        r   =   sqrt((xo(i)-xo(j))^2+(yo(i)-yo(j))^2);
        %% simple repulsive force node i if j is in the range of 
        if r < rpl_range,
            xn(i)= xn(i) + (xo(i) - xo(j));
            yn(i)= yn(i) + (yo(i) - yo(j));
            continue;
        end
        %% Attractivness based on Ligntness(or z) (z larger, lightness deemer)
        if  zo(i) > zo(j) && zo(j) < threshold,     % Brighter and smaller than 
            beta0   =   zo(j)/min_zo;               % Normalize Brightness
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

% TODO check if most of the point is moving in maximum step, if so, adjust
% the parameter
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
[xn,yn]=findrange(xn,yn,range);
end
%% Reduce the randomness during iterations
function alpha=newalpha(alpha,delta,threshold,zn)
    for i=1:numel(zn)
        if(zn(i) <= threshold)
            alpha(i)=alpha(i)*delta;
        else
            alpha(i) = alpha(i) / delta;
        end
        if(alpha(i) > 3)
            alpha(i) = 1;
        end
    end
end
% Make sure the fireflies are within the range
function [xn,yn]=findrange(xn,yn,range)

for i=1:length(yn),
   if xn(i)<=range(1,1), xn(i)=range(1,1); end
   if xn(i)>=range(1,2), xn(i)=range(1,2); end
   if yn(i)<=range(2,1), yn(i)=range(2,1); end
   if yn(i)>=range(2,2), yn(i)=range(2,2); end
end
%rn = xn.*xn + yn.*yn;
%for i = 1:length(yn)
%   if(rn(i) > 0.0004)
%       xn(i) = 0.02*xn(i) / (rn(i)^0.5) ;
%       yn(i) = 0.02*yn(i) / (rn(i)^0.5) ;
%   end
%end
end
%  ============== end =====================================
