clear all
close all
clc
%%

load('highres_numerical.mat'); % load numerical solution (from: Rotation_around_the_origin_center_unshifted.m)

rr=rr(:);
P_rot=P_rot(:);

%% Expression after rotation

ind1 = rr<= L0;
ind2 = rr>L0 & rr<R0;
ind3 = rr>=R0 & rr<=sqrt(R0^2+L0^2);

P_exact = zeros(size(rr));

r1 = rr(ind1);
r2 = rr(ind2);
r3 = rr(ind3);

P_exact(ind1) = 2.*R0./sqrt(R0.^2-r1.^2).*atan(sqrt( (R0-r1)./(R0+r1) ).* (sqrt(L0^2+R0^2)-R0)/L0) + ...
                R0./2./L0.*log( (L0.^2-(R0+r1).*(sqrt(L0^2+R0^2)-R0)) ./ (L0.^2-(R0-r1).*(sqrt(L0^2+R0^2)-R0)) )+...
                R0./sqrt(R0.^2-r1.^2).*(pi/2-atan(L0./sqrt(R0^2-r1.^2)));

P_exact(ind2) = 2.*R0./sqrt(R0.^2-r2.^2).*atan(sqrt( (R0-r2)./(R0+r2) ).* (sqrt(L0^2+R0^2)-R0)/L0) + ...
                R0./2./L0.*log( - ( (L0^2-(r2+sqrt(r2.^2-L0^2)).*(r2+R0)).*(L0^2-(r2+R0).*(sqrt(L0^2+R0^2)-R0)) ) ./ ( (L0^2-(r2+sqrt(r2.^2-L0^2)).*(r2-R0)).*(L0^2+(r2-R0).*(sqrt(L0^2+R0^2)-R0)) ) ) + ...
                R0./sqrt(R0^2-r2.^2).*( atan(sqrt((R0+r2)./(R0-r2)).*L0./(sqrt(r2.^2-L0^2)+r2)) + atan(sqrt((R0-r2)./(R0+r2)).*L0./(sqrt(r2.^2-L0^2)+r2)) - atan(sqrt((R0+r2)./(R0-r2))./L0.*(sqrt(R0.^2+L0^2)-R0)) - atan(sqrt((R0-r2)./(R0+r2))./L0.*(sqrt(R0.^2+L0^2)-R0)));          

P_exact(ind3) = R0./sqrt(r3.^2-R0.^2).*log(abs(R0./r3.*( sqrt(L0.^2+R0.^2)-R0+L0.*sqrt( (r3+R0)./(r3-R0) ) )./( sqrt(L0.^2+R0.^2)-R0-L0.*sqrt( (r3+R0)./(r3-R0) ) ) ) )+...
                R0./L0./2.*log( ( ( (r3+R0).*(r3+sqrt(r3.^2-L0^2)) - L0^2 ).* ( (r3+R0).*(sqrt(R0^2+L0^2) - R0) - L0^2 ) ) ./ ( (  (r3-R0).*(r3+sqrt(r3.^2-L0^2)) - L0^2 ).* ( (r3-R0).*(sqrt(R0^2+L0^2) - R0) + L0^2 ) ) )+...
                R0./2./sqrt(r3.^2-R0^2).*log(abs( (L0+sqrt(r3.^2-R0^2))./(-L0+sqrt(r3.^2-R0^2)).*(L0^3*R0-2.*L0.*R0.*r3.*(r3+sqrt(r3.^2-L0^2))+2.*r3.^2.*(r3+sqrt(r3.^2-L0^2)).*sqrt(r3.^2-R0.^2) - L0.^2.*(2.*r3+sqrt(r3.^2-L0^2)).*sqrt(r3.^2-R0^2))./(-L0^3*R0+2.*L0.*R0.*r3.*(r3+sqrt(r3.^2-L0^2))+2.*r3.^2.*(r3+sqrt(r3.^2-L0^2)).*sqrt(r3.^2-R0.^2) - L0.^2.*(2.*r3+sqrt(r3.^2-L0^2)).*sqrt(r3.^2-R0^2)) ));

P_exact=P_exact/R0/L0/pi;

figure, plot(rr,P_exact), hold on, plot(rr,P_rot),legend('exact','numerical');

figure, plot(rr,P_exact-P_rot), title('error')
