function [P_unique,rs_unique] = numerical_rotation(R0,L0,h,angle_step,M)
%numerical rotation of tent-like pdf.
if nargin <5
    M=200;
end
%% Placing a "tent" (Z-direction shift)
Res =  R0;

v0 = linspace(0,Res,M);
v = [-v0(end:-1:2) v0];
r = v;
z = v;

dr = mean(diff(r));
% dz = mean(diff(z));

[R,Z] = meshgrid(r,z);

enl = 0.5;
ind1 = abs(Z-h) <= L0./R0.*abs(R) & abs(R)<R0+enl*dr & abs(Z-h)<L0;
ind2 = abs(Z-h) > L0./R0.*abs(R)  & abs(R)<R0+enl*dr & abs(Z-h)<L0;

Pimg = zeros(size(ind1));

Pimg(ind1) = R0./(R0+abs(R(ind1))) ;
Pimg(ind2) = R0^2./(R0^2-abs(R(ind2)).^2).*(L0-abs(Z(ind2)-h))./L0;

Pimg=Pimg/2/R0/L0;

% sum(P(:))*dr*dz % should be approx 1

%% Rotating 
Pout = zeros(size(Pimg));
angles = 0:angle_step:360-angle_step;
for angle = angles
    Pout = Pout + imrotate(Pimg,angle,'bilinear','crop');
end
Pout = Pout/length(angles);

%% From radial image to 1d - data (r,P(r))

rs = sqrt(R(:).^2+Z(:).^2);
ind = rs>R0;
rs(ind) = [];
Pout(ind) =[];
tol = R0*1e-10;
[rs_unique,ia,~] = uniquetol(rs,tol);
P_unique = Pout(ia);

end

