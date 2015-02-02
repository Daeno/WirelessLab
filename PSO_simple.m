% ========================================================  
%                       PSO_simple                         
% ======================================================== 

function [best,M,countn,x,y,z] = PSO_simple(instr)
% n=number of fireflies
% MaxGeneration=number of pseudo time steps
if nargin<1,   instr=[20 300];     end
n=instr(1);  MaxGeneration=instr(2);

rand('seed',sum(100*clock));  % Reset the random generator
% ------ Four peak functions ---------------------
str1='exp(-(x-4)^2-(y-4)^2)+exp(-(x+4)^2-(y-4)^2)';
str2='+2*exp(-x^2-(y+4)^2)+2*exp(-x^2-y^2)+2*exp(-(x-1)^2-(y+2)^2)+7*exp(-(x-2)^2-y^2)';
funstr=strcat(str1,str2);
% Converting to an inline function
%f=vectorize(inline(funstr));

% range=[xmin xmax ymin ymax];
range=[-10 10 -10 10];
range2 = [-10 10;-10 10];
[x,y,z] = randfunc_g(range2,100);
% ------------------------------------------------
alpha  = 0.5;      % Randomness 0--1 (highly random)
gamma1 = 0.4;      % Absorption coefficient
gamma2 = 0.8;
delta  = 0.97;      % Randomness reduction (similar to 
                % an annealing schedule)
rigidity=0.8;
Vmax=2;
                
% ------------------------------------------------
% Grid values are used for display only
Ngrid=100;
dx=(range(2)-range(1))/Ngrid;
dy=(range(4)-range(3))/Ngrid;
[x,y]=meshgrid(range(1):dx:range(2),...
               range(3):dy:range(4));
%z=f(x,y);
% Display the shape of the objective function
figure(1);    surfc(x,y,z);


% ------------------------------------------------
% generating the initial locations of n fireflies
[xn,yn,vxn,vyn]=init_pso(n,range);
xb = xn;
yb = yn; 
lightb = interp2(x,y,z,xn,yn);
% Display the paths of fireflies in a figure with
% contours of the function to be optimized
 figure(2);
% Iterations or pseudo time marching
xg = 0; yg = 0; lightg = -30;

countn = 0;
%for i=1:MaxGeneration,     %%%%% start iterations
while(1)
    countn = countn + 1;
% Show the contours of the function
 contour(x,y,z,15); hold on;
% Evaluate new solutions
%zn=f(xn,yn);
zn = interp2(x,y,z,xn,yn);
xo = xn;
yo = yn;
% Ranking the fireflies by their light intensity
for i = 1:n,
    if zn(i) > lightb(i),
        lightb(i) = zn(i);
        xb(i) = xn(i);
        yb(i) = yn(i);
    end
end

[lightcurrent, index] = sort(lightb,'descend');
if lightcurrent(1) > lightg,
    lightg = lightcurrent(1);
    xg = xn(index(1)); yg = yn(index(1));    
end  
% Trace the paths of all roaming  fireflies
plot(xo,yo,'.','markersize',10,'markerfacecolor','g');
% Move all fireflies to the better locations
[xn,yn,vxn,vyn]=pso_move(xn,yn,vxn,vyn,xb,yb,xg,yg,alpha,gamma1,gamma2,range,Vmax,rigidity);
quiver(xo,yo,vxn,vyn,0);
drawnow;
M(countn) = getframe;

% Use "hold on" to show the paths of fireflies
    hold off;
    
% Reduce randomness as iterations proceed
alpha=newalpha(alpha,delta);
best = [xg,yg,lightg];
if(lightcurrent(numel(lightcurrent)) > 5) break; end;
end   %%%%% end of iterations
end



% ----- All subfunctions are listed here ---------
% The initial locations of n fireflies
function [xn,yn,vxn,vyn]=init_pso(n,range)
xrange=range(2)-range(1);
yrange=range(4)-range(3);
xn=rand(1,n)*xrange+range(1);
yn=rand(1,n)*yrange+range(3);
vxn=zeros(size(yn));
vyn=zeros(size(yn));
end

% Move all fireflies toward brighter ones
function [xn,yn,vxn,vyn]=pso_move(xn,yn,vxn,vyn,xb,yb,xg,yg,...
    alpha,gamma1,gamma2,range,Vmax,rigidity)
vxnn = rigidity*vxn+ (xb-xn)*gamma1*(1-rigidity) + (xg-xn)*gamma2*(1-rigidity);
vynn = rigidity*vyn+ (yb-yn)*gamma1*(1-rigidity) + (yg-yn)*gamma2*(1-rigidity);
v = vxnn.*vxnn + vynn.*vynn;
ni=size(yn,2);
for i=1:ni,
% The attractiveness parameter beta=exp(-gamma*r)
if v(i) > Vmax,
    vxnn(i)= 2*vxnn(i)/(v(i)^0.5);
    vynn(i)= 2*vynn(i)/(v(i)^0.5);
end
xn(i)=xn(i)+vxnn(i)+alpha.*(rand-0.5);
yn(i)=yn(i)+vynn(i)+alpha.*(rand-0.5);
end % end for i
vxn = vxnn;
vyn = vynn;
[xn,yn]=findrange(xn,yn,range);
xn
yn
xb
yb
xg
yg
end

% Reduce the randomness during iterations
function alpha=newalpha(alpha,delta)
alpha=alpha*delta;
end

% Make sure the fireflies are within the range
function [xn,yn]=findrange(xn,yn,range)
for i=1:length(yn),
   if xn(i)<=range(1), xn(i)=range(1); end
   if xn(i)>=range(2), xn(i)=range(2); end
   if yn(i)<=range(3), yn(i)=range(3); end
   if yn(i)>=range(4), yn(i)=range(4); end
end
end
%  ============== end =====================================


