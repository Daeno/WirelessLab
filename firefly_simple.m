% ======================================================== % 
% Files of the Matlab programs included in the book:       %
% Xin-She Yang, Nature-Inspired Metaheuristic Algorithms,  %
% Second Edition, Luniver Press, (2010).   www.luniver.com %
% ======================================================== %    


% =========================================================% 
% Firefly Algorithm by X S Yang (Cambridge University)     %
% Usage: firefly_simple([number_of_fireflies,MaxGeneration])
%  eg:   firefly_simple([12,50]);                          %
% ======================================================== %
% This is a demo for 2D functions; for higher dimenions,   %
% you should use fa_ndim.m or fa_mincon.m                  %
% Parameters choice: 
% Gamma should be linked with scales. Otherwise, the FA    %
% the efficiency will be significantly reduced because     %
% the beta term may be too small.                          %
% Similarly, alpha should also be linked with scales,      %
% the steps should not too large or too small, often       %
% steps are about 1/10 to 1/100 of the domain size.        %
% In addition, alpha should be reduced gradually           %
% using alpha=alpha_0 delta^t during eteration t.          %
% Typically, delta=0.9 to 0.99 will be a good choice.      %
% ======================================================== %

function [N, count_n,x,y,z] = firefly_simple(instr,x,y,z)
% n=number of fireflies
% MaxGeneration=number of pseudo time steps
if nargin<1,   instr=[100 10000];     end
n=instr(1);  MaxGeneration=instr(2);
% Show info
help firefly_simple.m
%rand('state',0);  % Reset the random generator
% ------ Four peak functions ---------------------
str1='3*exp(-(x-5)^2-(y-4)^2) - 4*exp(-(x+2)^2-(y-1.5)^2)';
str2='+6*exp(-x^2-(y+5)^2) - 0.2*exp(-x^2-y^2)';
str3 = '- 0.3 * cos(0.7*3.14159*x) - 0.4*cos(0.35*3.14159*y) - 0.3 * cos(3*3.14159*x) - 0.4*cos(4*3.14159*y)+0.7'
funstr=strcat(str1,str2);
funstr = strcat(funstr,str3);
% Converting to an inline function
f=vectorize(inline(funstr));
% range=[xmin xmax ymin ymax];
range=[-10 10 ; -10 10];

% ------------------------------------------------
%alpha=0.2;      % Randomness 0--1 (highly random)
alpha = 0.2 * ones(1,instr(1));
gamma=1;      % Absorption coefficient
delta=0.7;      % Randomness reduction (similar to 
                % an annealing schedule)
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
N = 0;
% ------------------------------------------------
% generating the initial locations of n fireflies
[xn,yn,Lightn]=init_ffa(n,range);
% Display the paths of fireflies in a figure with
% contours of the function to be optimized
 figure(2);
 count_n = 0;
% Iterations or pseudo time marching

while(count_n < MaxGeneration) %for i=1:MaxGeneration,     %%%%% start iterations
% Show the contours of the function
%axis equal
count_n = count_n + 1;
%contour(x,y,z,15); hold on;
% Evaluate new solutions
%zn=f(xn,yn);
zn = interp2(x,y,z,xn,yn);
% Ranking the fireflies by their light intensity
[Lightn,Index]=sort(zn);

N(count_n) = sum(Lightn > 2.5);

xn=xn(Index); yn=yn(Index); alpha = alpha(Index);
    
xo=xn;   yo=yn;    Lighto=Lightn;
% Trace the paths of all roaming  fireflies
%plot(xn,yn,'.','markersize',10,'markerfacecolor','g');

% Move all fireflies to the better locations
alpha=newalpha(alpha,delta,zn);
[xn,yn]=ffa_move(xn,yn,Lightn,xo,yo,Lighto,alpha,gamma,range);
%quiver(xo,yo,xn-xo,yn-yo);
%drawnow;
%M(count_n) = getframe;
% Use "hold on" to show the paths of fireflies
%    hold off;
    
% Reduce randomness as iterations proceed
if(Lightn(1) > 2.5) break;  end
end   %%%%% end of iterations


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
    
    %if Lightn(i) > 2.5, 
            %continue; 
    %end
% The attractiveness parameter beta=exp(-gamma*r)
    for j=1:nj,
        r=sqrt((xn(i)-xo(j))^2+(yn(i)-yo(j))^2);
        if r < 1,
                xn(i)= xn(i)+(xn(i)-xo(j))+alpha(i).*(rand-0.5)*0.25;
                yn(i)= yn(i)+(yn(i)-yo(j))+alpha(i).*(rand-0.5)*0.25;
                continue;
        end
        if Lightn(i)<Lighto(j) && Lighto(j) > 2.5, % Brighter and more attractive
            beta0=Lightn(j)/Lightn(length(Lightn));     beta=beta0*exp(-gamma*r.^2);
            %*
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
            alpha(i) = alpha(i) / delta;
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
%rn = xn.*xn + yn.*yn;
%for i = 1:length(yn)
%   if(rn(i) > 0.0004)
%       xn(i) = 0.02*xn(i) / (rn(i)^0.5) ;
%       yn(i) = 0.02*yn(i) / (rn(i)^0.5) ;
%   end
%end
end
%  ============== end =====================================
