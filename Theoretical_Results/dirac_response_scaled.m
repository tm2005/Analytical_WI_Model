function [y] = dirac_response_scaled(r,h,L0,R0)
% Rotation of dirac line

y = real( 1./(sqrt(r.^2-h.^2)) );

% scaling
y = y*real(asin(R0/h))/2/R0/L0;
