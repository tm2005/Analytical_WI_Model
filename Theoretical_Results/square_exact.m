function [y] = square_exact(r,h,L0,R0)
% solution when rotation a square approx of pdf
y = real( asin( (h+L0)./r) - asin( (h-L0)./r));

if h == L0
      y = real( asin( (h+L0)./r ) );
end
y = y/4/L0/R0/pi;
end