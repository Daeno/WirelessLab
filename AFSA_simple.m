function [best,count_n,M,x,y,z] = AFSA_simple(instr,x,y,z)
    % n=number of aitifical fishes
    if nargin<1,   instr=[60 50];     end
    n=instr(1);  MaxGeneration=instr(2);
    % Show info
    range=[-10 10 ;-10 10];

    % ------------------------------------------------
    eta = 0.01;      % stagnation criteria
    delta=0.3;      % visual coefficient
    % ------------------------------------------------
    % Grid values are used for display only
    %Ngrid=100;
    %dx=(range(2)-range(1))/Ngrid;
    %dy=(range(4)-range(3))/Ngrid;
    %[x,y]=meshgrid(range(1):dx:range(2),...
    %               range(3):dy:range(4));
    %z=f(x,y);
    % Display the shape of the objective function
    %figure(1);    surfc(x,y,z);
    if(nargin < 4)
        [x,y,z] = randfunc_g(range,100);
    end
    figure(1);
    surfc(x,y,z);

    % ------------------------------------------------
    % generating the initial locations of n fireflies
    [xn,yn,foodn]=init_afsa(n,range);
    % Display the paths of fireflies in a figure with
    % contours of the function to be optimized
    figure(2);
    count_n = 0;
    % Iterations or pseudo time marching
    prev_food = zeros(size(xn));
    sat_cnt = zeros(size(xn));
    while(1) %for i=1:MaxGeneration,     %%%%% start iterations
        % Show the contours of the function
        axis equal
        count_n = count_n + 1;
        contour(x,y,z,15); hold on;
        % Evaluate new solutions
        %zn=f(xn,yn);
        zn = interp2(x,y,z,xn,yn);
        % Ranking the fireflies by their light intensity
        [foodn,Index]=sort(zn);
        xn=xn(Index); yn=yn(Index);
    
        xo=xn;   yo=yn;    foodo=foodn;
        % Trace the paths of all roaming  fireflies
        plot(xn,yn,'.','markersize',10,'markerfacecolor','g');
        
        % Move all fireflies to the better locations
        %alpha=newalpha(alpha,delta,zn);
        [xn,yn,prev_food,sat_cnt]=afsa_move(xn,yn,xo,yo,foodo,prev_food,sat_cnt,eta,delta,range,x,y,z);
        quiver(xo,yo,xn-xo,yn-yo,0);
        drawnow;
        M(count_n) = getframe;
        % Use "hold on" to show the paths of fireflies
        hold off;
        if(foodo(1) > 2.5) break;  end
        % Reduce randomness as iterations proceed
    
    end   %%%%% end of iterations
    best(:,1)=xo'; best(:,2)=yo'; best(:,3)=foodo';

end
% ----- All subfunctions are listed here ---------
% The initial locations of n fireflies
function [xn,yn,foodn]=init_afsa(n,range)
    xrange=range(1,2)-range(1,1);
    yrange=range(2,2)-range(2,1);
    xn=rand(1,n)*xrange+range(1,1);
    yn=rand(1,n)*yrange+range(2,1);
    foodn=zeros(size(yn));
end

% Move fishes
function [xn,yn,prev_food, sat_cnt]=afsa_move(xn,yn,xo,yo,foodo,prev_food,sat_cnt,eta,delta,range,x,y,z)
    visual = 2;%delta * max(range(:,2) - range(:,1));
    [nums,nearby, cp, best_food, best_i] = count_visual(xo,yo,foodo,visual);
    %visual
    %nums
    %nearby
    %cp
    %best_food
    %best_i
    nxn = xn;
    nyn = yn;
    for i = 1:size(xn,2)
        if(nums(i) == 0)
            %fprintf('no fish\n');
            [nxn(i), nyn(i)] = afsa_Rnd(xo(i),yo(i),range,0.02);
        else
            if(nums(i) > size(xn,2) * delta)
                %fprintf('crowed\n');
                % search
                %[nxn(i), nyn(i)] = afsa_Srh(i,xo,yo,foodo,nums(i),nearby(:,i),visual,range,x,y,z);
                [nxn(i), nyn(i)] = afsa_Srh(i,xo,yo,foodo,nums,nearby(:,i),range,visual,x,y,z);
            else
                if(interp2(x,y,z,cp(1,i),cp(2,i)) > foodo(i))
                    %fprintf('swarm\n');
                    [nxn(i), nyn(i)] = afsa_move_toward(xo(i),yo(i),cp(1,i),cp(2,i));
                else
                    %fprintf('chase or search\n');
                    tempx = zeros(2,1);
                    tempy = zeros(2,1);
                    %%% central point preprocessing
                    for i=1:size(xo,2)
                        cp(:,i) = cp(:,i) / nums(i);
                    end
                    eva1 = interp2(x,y,z,cp(1,i),cp(2,i));
                    if(eva1 > foodo(i))
                        %fprintf('      swarm,');
                        tempx(1,1) = cp(1,i);
                        tempy(1,1) = cp(2,i);
                    else
                        [tempx(1,1), tempy(1,1)] = afsa_Srh(i,xo,yo,foodo,nums,nearby(:,i),range,visual,x,y,z);
                    end
                    eva2 = 0;
                    dx = 0;
                    dy = 0;
                    for iter = 1:5
                        dx = rand-0.5;
                        dy = rand-0.5;
                        l = rand;
                        d = (dx*dx + dy*dy)^0.5;
                        dx = dx*l/d*visual;
                        dy = dy*l/d*visual;
                        eva2 = interp2(x,y,z,xo(i)+dx,yo(i)+dy);
                        if(eva2 > foodo(i))
                            %fprintf('   chase,');
                            tempx(2,1) = xo(i) + dx;
                            tempy(2,1) = yo(i) + dy;
                        else
                            [tempx(2,1), tempy(2,1)] = afsa_Srh(i,xo,yo,foodo,nums,nearby(:,i),range,visual,x,y,z);
                        end
                    end
                    if(eva2 > eva1)
                        nxn(i) = tempx(2,1);
                        nyn(i) = tempy(2,1);
                    else
                        nxn(i) = tempx(1,1);
                        nyn(i) = tempy(1,1);
                    end
                end
            end
        end
    end
    for i = 1:size(xn,2)
        tmp_food = interp2(x,y,z,nxn(i),nyn(i));
        if(foodo(i) > tmp_food)
            xn(i) = xo(i);
            yn(i) = yo(i);
            sat_cnt(i) = sat_cnt(i) + 1;
        else
            xn(i) = nxn(i);
            yn(i) = nyn(i);
            if(abs(prev_food(i) - tmp_food) / tmp_food < eta)
                sat_cnt(i) = sat_cnt(i) + 1;
            else
                sat_cnt(i) = 0;
            end
        end
        
        if(sat_cnt(i) == 3)
            sat_cnt(i) = 2;
            [xn(i), yn(i)] = afsa_Rnd(xo(i),yo(i),range,0.02);
        end
        prev_food(i) = interp2(x,y,z,xn(i),yn(i));
    end
end

function [nums,nearby, cp, best_food, best_i]= count_visual(xo,yo,foodo,visual)
    m = size(xo,2);
    nums = zeros(size(xo));
    nearby = zeros(size(xo));
    cp = zeros(2,size(xo,2));
    %cp(1,:) = xo;
    %cp(2,:) = yo;
    best_food = nums;
    best_i = nums;
    visual_sq = visual*visual;
    for i = 1:m
        for j = 1:m
            if((xo(j) - xo(i))^2 + (yo(j) - yo(i))^2 < visual_sq && i ~= j)
                nums(i) = nums(i) + 1;
                if(best_food(i) < foodo(j))
                    best_food(i) = foodo(j);
                    best_i(i) = j;
                end
                nearby(nums(i),i) = j;
                cp(1,i) = cp(1,i) + xo(j);
                cp(2,i) = cp(2,i) + yo(j);
            end
        end
    end
end

function [x, y] = afsa_Rnd(xo_,yo_,range,visual)
    l1 = rand;
    l2 = rand;
    if(l1 > 0.5)
        if(range(1,2) - xo_ > visual)
            x = xo_ + l2*visual;
        else
            x = xo_ + l2*(range(1,2) - xo_);
        end
    else
        if(xo_ - range(1,1) > visual)
            x = xo_ - l2 * visual;
        else
            x = xo_ - l2 * (xo_ - range(1,1));
        end
    end
    
    l1 = rand;
    l2 = rand;
    if(l1 > 0.5)
        if(range(1,2) - yo_ > visual)
            y = yo_ + l2*visual;
        else
            y = yo_ + l2*(range(1,2) - yo_);
        end
    else
        if(yo_ - range(1,1) > visual)
            y = yo_ - l2 * visual;
        else
            y = yo_ - l2 * (yo_ - range(1,1));
        end
    end
end

function [xn, yn] = afsa_Srh(idx,xo,yo,foodo,nums,nearby_,range,visual,x,y,z)
    best = foodo(idx);
    best_i = idx;
    perm = randperm(nums(idx));
    for i = 1:nums(idx)
        if(best < foodo(nearby_(perm(i))))
            best = foodo(nearby_(perm(i)));
            best_i = nearby_(perm(i));
        end
    end
    if(best_i == idx)
        for i = 1:5
            [xn, yn] = afsa_Rnd(xo(idx),yo(idx),range,visual);
            if(interp2(x,y,z,xn,yn) > foodo(idx)) break; end
        end
    else
%function [x, y] = afsa_Srh(idx,xo,yo,nearby)
%    id = ceil(rand * size(nearby,1));
%    id = nearby(id);
        [xn, yn] = afsa_move_toward(xo(idx),yo(idx),xo(best_i),yo(best_i));
    end
end

function [x, y] = afsa_move_toward(xo_, yo_, xn_, yn_)
    l = rand;
    if((xn_ - xo_)^2 + (yn_ - yo_)^2 > 0.0004)
        dx = xn_ - xo_;
        dy = yn_ - yo_;
        x = xo_ + 0.02*l*dx/((dx^2 + dy^2)^0.5);
        y = yo_ + 0.02*l*dy/((dx^2 + dy^2)^0.5);
    else
        x = xo_ + l*(xn_ - xo_);
        y = yo_ + l*(yn_ - yo_);
    end
end
%  ============== end =====================================
