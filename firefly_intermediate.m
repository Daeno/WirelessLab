function [ N,count_n,x,y,z ] = firefly_intermediate( instr,x,y,z )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin<1,   instr=[200 50];     end
n=instr(1);  MaxGeneration=instr(2);
best = zeros(1,instr(2));
% Show info
help firefly_simple.m
%rand('seed',sum(100*clock)); % Reset the random generator

% range=[xmin xmax ymin ymax];
range=[-10 10 ; -10 10];

%%%%%%%
%
%crowded_pattern=[-0.1 -0.2 -0.2 -0.1;
%                 -0.2  -0.5 -0.5 -0.2;
%                 -0.2  -0.5 -0.5 -0.2;
%                 -0.1 -0.2 -0.2 -0.1];
%             
str='-0.25*exp(0.05*(-x^2-y^2))';
f=vectorize(inline(str)); 
[xt,yt]=meshgrid(-3:0.5:3,...
               -3:0.5:3);
crowded_pattern_exp=f(xt,yt);
             
%%%%%%%             

% ------------------------------------------------
%alpha=0.2;      % Randomness 0--1 (highly random)
alpha = 0.2 * ones(1,instr(1));
gamma=1;      % Absorption coefficient
delta=0.7;      % Randomness reduction (similar to 
                % an annealing schedule)
% ------------------------------------------------


% Display the shape of the objective function
%figure(1);    surfc(x,y,z);
if(nargin < 4)
    [x,y,z] = randfunc_g(range,100);
end
%z = zeros(size(z));
figure(1);
surfc(x,y,z);

% ------------------------------------------------
% generating the initial locations of n fireflies
[xn,yn,Lightn]=init_ffa(n,range);
% Display the paths of fireflies in a figure with
% contours of the function to be optimized
 figure(2);
 count_n = 0;
% Iterations or pseudo time marching

for i=1:MaxGeneration,     %%%%% start iterations
    % Show the contours of the function
    axis equal
    count_n = count_n + 1;

%%%%%%%%    
    N(count_n) = sum(Lightn > 2.5);
    ztemp = z;
    % Evaluate new solutions
    %zn=f(xn,yn);
    for i = 1:n,
        rx = 1;
        ry = 1;
        for j = 1:size(x,2)
            for k = 1:size(y,1)
            if(x(k,j) - xn(i) > 0) continue; 
            else
                rx = j;
            end
            if(y(k,j) - yn(i) > 0) continue;
            else
                ry = k;
            end
            end
        end
        ry;
        if (rx < 1 || ry < 1) continue; end
        for p = rx-6:rx+6,
            if (p > size(x,1) || p < 1), continue; end
            for q = ry-6:ry+6,
                if (q > size(x,1) || q < 1), continue; end
                ztemp(q,p)= ztemp(q,p) + crowded_pattern_exp(q - ry + 7,p-rx+7);
            end
        end
    end
    %figure(1);
    %surfc(x,y,ztemp);
    %figure(2);
    %contour(x,y,ztemp,15); hold on;
    
%%%%%%%%
    
    zn = interp2(x,y,ztemp,xn,yn);
    % Ranking the fireflies by their light intensity
    [Lightn,Index]=sort(zn);
    xn=xn(Index); yn=yn(Index); alpha = alpha(Index);
    xo=xn;   yo=yn;    Lighto=Lightn;
    % Trace the paths of all roaming  fireflies
    plot(xn,yn,'.','markersize',10,'markerfacecolor','g');
    
    
    
    % Move all fireflies to the better locations
    alpha=newalpha(alpha,delta,zn);
    [xn,yn]=ffa_move(xn,yn,Lightn,xo,yo,Lighto,alpha,gamma,range);
    %quiver(xo,yo,xn-xo,yn-yo);
    %drawnow;
    %M(count_n) = getframe;
    if(Lightn(1) > 2.5) break;  end
    % Use "hold on" to show the paths of fireflies
    hold off;
    
    
end   %%%%% end of iterations
%best(:,1)=xo'; best(:,2)=yo'; best(:,3)=Lighto';

end
% ----- All subfunctions are listed here ---------
% The initial locations of n fireflies
function [xn,yn,Lightn]=init_ffa(n,range)
    xrange=range(1,2)-range(1,1);
    yrange=range(2,2)-range(2,1);
    xn=rand(1,n)*xrange+range(1,1);
    yn=rand(1,n)*yrange+range(2,1);
    Lightn=zeros(size(yn));
end

% Move all fireflies toward brighter ones
function [xn,yn]=ffa_move(xn,yn,Lightn,xo,yo,...
    Lighto,alpha,gamma,range)
    ni=size(yn,2); nj=size(yo,2);
    for i=1:ni,
    % The attractiveness parameter beta=exp(-gamma*r)
   
%%%%%%%    
    
        if Lightn(i) > 2.5, 
            continue; 
        end
        
%%%%%%%        
        
        for j=1:nj,
            r=sqrt((xn(i)-xo(j))^2+(yn(i)-yo(j))^2);
            if Lightn(i)<Lighto(j), % Brighter and more attractive
                beta0=Lightn(j)/Lightn(length(Lightn));     beta=beta0*exp(-gamma*r.^2);
                xn(i)=xn(i).*(1-beta)+xo(j).*beta+alpha(i).*(rand-0.5)*0.25;
                yn(i)=yn(i).*(1-beta)+yo(j).*beta+alpha(i).*(rand-0.5)*0.25;
            end
        end % end for j
    end % end for i
   
dx = xn - xo;
dy = yn - yo;
r = dx.*dx + dy.*dy;
for idx = 1:numel(dx)
   if(r(idx) > 0.0004)
        dx(idx) = dx(idx) / r(idx).^0.5 * 0.2;
        dy(idx) = dy(idx) / r(idx).^0.5 *0.2;
   end
end
xn = dx + xo;
yn = dy + yo;
    
    [xn,yn]=findrange(xn,yn,range);
end
% Reduce the randomness during iterations
function alpha=newalpha(alpha,delta,zn)
    for i=1:numel(zn)
        if(zn(i) >= 2)
            alpha(i)=alpha(i)*delta;
        else
            alpha(i)=alpha(i)/delta;
        end
        if(alpha(i) > 1)
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
end
%  ============== end =====================================


