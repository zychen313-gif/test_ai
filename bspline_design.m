function B = bspline_design(x, xmin, xmax, order, nKnots)
x = max(min(x, xmax), xmin);
kn = linspace(xmin, xmax, nKnots);
t = [repmat(xmin,1,order), kn(2:end-1), repmat(xmax,1,order)];
nb = numel(t) - order - 1;
B = zeros(numel(x), nb);
for i = 1:nb
    B(:,i) = cox_de_boor(x, i, order, t);
end
B = B ./ max(sum(B,2), eps);
end

function y = cox_de_boor(x, i, k, t)
if k == 1
    y = double(x >= t(i) & x < t(i+1));
    if i == numel(t)-1
        y(x == t(end)) = 1;
    end
else
    d1 = t(i+k-1) - t(i);
    d2 = t(i+k) - t(i+1);
    y1 = zeros(size(x)); y2 = zeros(size(x));
    if d1 > 0
        y1 = (x - t(i)) ./ d1 .* cox_de_boor(x, i, k-1, t);
    end
    if d2 > 0
        y2 = (t(i+k) - x) ./ d2 .* cox_de_boor(x, i+1, k-1, t);
    end
    y = y1 + y2;
end
end
