function zPlane = three_point_level_plane(x, y, z)
r0 = 75;
angles = [0, 120, 240];
px = r0 * cosd(angles);
py = r0 * sind(angles);
F = scatteredInterpolant(x, y, z, 'natural', 'nearest');
pz = F(px(:), py(:));
A = [px(:), py(:), ones(3,1)];
coef = A \ pz(:);
zPlane = coef(1)*x + coef(2)*y + coef(3);
end
