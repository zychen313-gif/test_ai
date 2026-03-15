function d = load_and_preprocess_wafer(filePath, cfg)
raw = fun_read_wafer_map(filePath);
x = raw(:,6); y = raw(:,7); z = raw(:,8);
valid = isfinite(x) & isfinite(y) & isfinite(z);
x = x(valid); y = y(valid); z = z(valid);

mu = mean(z); sd = std(z);
outlier = abs(z - mu) > cfg.outlierSigma * max(sd, eps);
if any(outlier) && nnz(~outlier) >= 20
    F = scatteredInterpolant(x(~outlier), y(~outlier), z(~outlier), 'natural', 'nearest');
    z(outlier) = F(x(outlier), y(outlier));
end

r = hypot(x, y);
th = mod(atan2d(y, x), 360);
notchMask = angular_distance_deg(th, cfg.notchCenterDeg) <= cfg.notchHalfWidthDeg;
edgeMask = (r > cfg.rMin) & (r <= cfg.rMax) & ~notchMask;

zPlane = three_point_level_plane(x, y, z);
zY = z - zPlane;

d.filePath = filePath;
d.x = x; d.y = y; d.r = r; d.theta = th;
d.zX = z;
d.zY = zY;
d.edgeMask = edgeMask;
end
