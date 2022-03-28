function [x,y] = nema_nu4_five(Nevents,R0,phi0,r0,xshift,yshift)
% NEMA Phantom with small active cricles:
% Nevents (int) - total number of events
% R0 - Distance from (xshift, yshift) to centers of active circles
% phi0 - angle between x-axis and the smallest active circle
% r0 (vector) - radius of all active circles

    N0 = round(1.7*Nevents);
    Ncircles = length(r0);
    phistep = 2*pi/Ncircles;
    r0_norm = sum(r0.^2);
    x=[];
    y=[];
    for i = 1:Ncircles
        N = round(N0*r0(i)^2/r0_norm);
%         if r0(i) <=2
%             N=N+N*
        xc = xshift + R0*cos(phi0 + phistep*(i-1));
        yc = yshift + R0*sin(phi0 + phistep*(i-1));
        
        xt = (2*r0(i))*rand(N,1) - r0(i) +xc;
        yt = (2*r0(i))*rand(N,1) - r0(i) +yc;
        
        ind = (xt  - xc).^2 + (yt - yc).^2 >= r0(i)^2;
        
        xt(ind)=[];
        yt(ind)=[];
        x = [x;xt];
        y = [y;yt];
    
    end
    ind = randperm(length(x),length(x)-Nevents);
    x(ind)=[];
    y(ind)=[];

end

