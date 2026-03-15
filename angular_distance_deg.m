function d = angular_distance_deg(a, b)
d = abs(mod(a - b + 180, 360) - 180);
end
