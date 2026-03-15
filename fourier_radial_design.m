function B = fourier_radial_design(r, tDeg)
rr = (r - 100) / 50;
t = deg2rad(tDeg);
B = [ones(size(rr)), rr, rr.^2, rr.^3, ...
     cos(t), sin(t), cos(2*t), sin(2*t), cos(3*t), sin(3*t), ...
     rr.*cos(t), rr.*sin(t), rr.*cos(2*t), rr.*sin(2*t)];
end
