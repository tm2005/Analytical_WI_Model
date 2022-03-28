function [I,R] = ideal_shifted(N,R0t,R0,L0 ,h)
% Not approximatied PDF 
    I =  zeros(N,N);
    x = linspace(-R0t,R0t,N);

    [R,Z] = meshgrid(x,x);

    ind1 = abs(Z-h) <= L0./R0.*abs(R) & abs(R)<R0 & abs(Z-h)<L0;% & abs(R)<=R0 & abs(Z)<=L0;
    ind2 = abs(Z-h) > L0./R0.*abs(R) & abs(R)<R0 & abs(Z-h)<L0; %& abs(Z)<=L0 & abs(R)<=R0;
   
    
    I(ind1) = R0./(R0+abs(R(ind1))) ;
    I(ind2) = R0^2./(R0^2-abs(R(ind2)).^2).*(L0-abs(Z(ind2)-h))./L0;


end