function [y] = triangle_exact(r,h,L0,R0)
% solution when rotation a trinagular approx of pdf
if h~=0
      y = real( (L0+h).*asin((L0+h)./r) - 2.*h.*asin(h./r) + (L0-h).*asin((L0-h)./r) +...
        sqrt(r.^2-(L0+h).^2) - 2.*sqrt(r.^2-h.^2) + sqrt(r.^2-(L0-h).^2) );
end

if h==L0
      y = real( (L0+h).*asin((L0+h)./r) - 2.*h.*asin(h./r) + 0 +...
        sqrt(r.^2-(L0+h).^2) - 2.*sqrt(r.^2-h.^2) + sqrt(r.^2-(L0-h).^2) );
end

if h == 0
      y = real( 2*(L0).*asin((L0)./r) + 2*sqrt(r.^2-(L0).^2) - 2.*r);
end
y = y/2/L0^2/R0/pi;
end