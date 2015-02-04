function [N, count_n, x, y, z] = baseline_re_search(instr,x,y,z, p_x, p_y)
%% Parameter declaration
if nargin<1,   instr=[200 50 1];     end
n=instr(1);  MaxGeneration=instr(2); draw = instr(3);

range=[-10 10 ; -10 10];    % range=[xmin xmax ymin ymax];
N = 0;

step = 0.2;                 % maximum single step
threshold = -1;
% Display the shape of the objective function

if(nargin < 4)
[x,y,z] = randfunc_g(range,100);
end
figure(1);
surfc(x,y,z);

% generating the initial locations of n fireflies
if(nargin < 6)
    [xn, yn]=init_pos(n,range(1,:),range(2,:));
else
    xn = p_x;
    yn = p_y;
end

% constants for reduction of intensity
[m, n] = meshgrid(-0.4:0.2:0.4, -0.4:0.2:0.4);
r = m.^2 + n.^2;
a = 2.*power(4, -r);

% compensate discoveries
cntd = zeros(size(xn));
numd = zeros(size(xn));
mapd = zeros(n, size(x));


figure(2);
count_n = 0;                    % Loop Counter

%% Main Loop
while(count_n < MaxGeneration)
    count_n = count_n + 1;
    xo = xn;
    yo = yn;
    zo = interp2(x,y,z,xo,yo);              % Evaluate new solutions
    [Lightn, Index]=sort(zo);               % Ranking the fireflies by their light intensity

    N(count_n) = sum(Lightn < threshold);

    [xn,yn]= simple_move(x,y,z,xo,yo,zo,step,range, threshold); % Move all fireflies to the better locations
    if(draw)
        axis equal;
        contour(x,y,z,15);
        hold on;
        plot(xo(zo>threshold),yo(zo>threshold),'.','markersize',10,'markerfacecolor','r');
        plot(xo(zo<=threshold),yo(zo<=threshold),'g.','markersize',10,'markerfacecolor','r');
        quiver(xo,yo,xn-xo,yn-yo);
        drawnow;
        hold off;
    end
    
    if(Lightn(n) < threshold) break;  end
end % Main Loop End

end
%% ----- All subfunctions are listed here ---------
% subfunction used in main function
%% Move all fireflies toward brighter ones
function [xn,yn]= simple_move(x,y,z,xo,yo,zo,step,range, threshold)
    ni=size(yo,2);
    xn = xo;
    yn = yo;
    zn = zo;
    for i=1:ni,
        %% Stop when small enough
        if zo(i) < threshold, 
            %continue;
            % count to 3 and set flag to discovered
            if cntd(i) < 3,
                cntd(i) = cntd(i) + 1;
            else
                cntd(i) = 0;
                % first colume denotes the numbered of flag discovered
                % the rest of colume records the flag position.
                numd(i) = numd(i) + 1;
                [I, J] = findindex(xn, yn, range, step);
                gridrange = size(x);
                for m = I-3:I+3,
                    for n = J-3:J+3,
                        if m > 0 && m < gridrange(1, 1) && n > 0 && n < gridrange(1, 2),
                            mapd(i,m,n) = mapd(i,m,n) + a(m-I+4, n-J+4);
                        end
                    end
                end
            end
        end
        %% generate next position
        r = rand * step;
        theta = rand * 2 * pi;
        xn(i) = xo(i) + r*cos(theta);
        yn(i) = yo(i) + r*sin(theta);
        [xn(i),yn(i)]=findrange(xn(i),yn(i),range);
        zn(i) = interp2(x,y,z,xn(i),yn(i));
        %% if new z is bigger than original z, pick a new one
        if (zn(i) > zo(i))                   
            xn(i) = xo(i) - r*cos(theta);
            yn(i) = yo(i) - r*sin(theta);
            [xn(i),yn(i)]=findrange(xn(i),yn(i),range);
        end
    end % end for i
    
end

function [I, J] = findindex(p_x, p_y, range, step)
    I = ceil((p_x-range(1, 1))/step);
    J = ceil((p_y-range(2, 1))/step);
end


