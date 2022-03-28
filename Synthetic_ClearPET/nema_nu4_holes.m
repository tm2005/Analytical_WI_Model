function [x,y] = nema_nu4_holes(Nevents,R0,phi0,Rc,r0,xshift,yshift)
% NEMA Phantom with holes:
% Nevents (int) - total number of events
% R0 - Distance from (xshift, yshift) to centers of active circles
% phi0 - angle between x-axis and the smallest active circle
% r0 (vector) - radius of all active circles

    Ncircles = length(r0);
    phistep = 2*pi/Ncircles;
    factor0 = sum(r0.^2)/R0^2;
    N0 = round(1.5*4/(pi*(1-factor0))*Nevents);
    x = (2*R0)*rand(N0,1) - R0 +xshift;
    y = (2*R0)*rand(N0,1) - R0 +yshift;
    ind = (x  - xshift).^2 + (y - yshift).^2 >= R0^2;
    x(ind)=[];
    y(ind)=[];    
    
    for i = 1:Ncircles
        xc = xshift + Rc*cos(phi0 + phistep*(i-1));
        yc = yshift + Rc*sin(phi0 + phistep*(i-1));
        
        ind = (x  - xc).^2 + (y - yc).^2 <= r0(i)^2;
        
        x(ind)=[];
        y(ind)=[];
    
    end
    x(Nevents+1:end)=[];
    y(Nevents+1:end)=[];

end
