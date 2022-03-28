function [P] = rot_trian_1d_Rij(r,hij,Rij,Lij)
% Expression for one component (one crystal-to-crystal response) using triangular approx
R0s=Rij;
L0s=Lij;
h=hij;
C = 1./2./L0s./R0s./pi./L0s;
C=C.*L0s^2; % Appendix D
indr0 = r==0;
P = C.*real( (L0s+h).*asin((L0s+h)./r) - 2*h*asin(h./r) + (L0s-h).*asin((L0s-h)./r) + sqrt(r.^2-(L0s+h).^2) - 2.*sqrt(r.^2-h.^2) + sqrt(r.^2-(L0s-h).^2) );
if h==0
    P = C.*real( (L0s).*asin((L0s)./r) + (L0s).*asin((L0s)./r) + sqrt(r.^2-(L0s).^2) - 2.*r + sqrt(r.^2-(L0s).^2) );
    P(indr0) = C.*real( (L0s+h).*pi/2 - 2*h*pi/2 + (L0s-h).*pi/2+ sqrt(0.^2-(L0s+h).^2) - 2.*sqrt(0.^2-h.^2) + sqrt(0.^2-(L0s-h).^2) );
end
end
