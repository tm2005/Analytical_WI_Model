function y = rescale01(x)
y = zeros(size(x));
xm = min(x(:));
xM = max(x(:));
y = (x - xm)/(xM-xm);
end