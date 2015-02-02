% ­°§C©Y«×
function [Xq Yq Vq] = randfunc_g(range,nstep)
    
    %range = [-6 6 -6 6];
    [X Y] = meshgrid(range(1,1):2:range(1,2),range(2,1):2:range(2,2));
    V = rand(size(X));
    for i=1:20
        V(ceil(rand(1) * numel(V))) = 9;
    end
    [Xq Yq] = meshgrid(range(1,1):(range(1,2) - range(1,1))/nstep:range(1,2),...
        range(2,1):(range(2,2) - range(2,1))/nstep:range(2,2)) ;
    Vq = interp2(X,Y,V,Xq,Yq,'spline');
    
    %figure; 
    %surfc(Xq,Yq,Vq);
end