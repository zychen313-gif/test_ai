function B = zernike_design(r, thetaDeg, rNormMax, nMax)
rho = min(max(r ./ rNormMax, 0), 1);
theta = deg2rad(thetaDeg);
terms = zernike_term_list(nMax);
B = zeros(numel(r), size(terms,1));
for k = 1:size(terms,1)
    n = terms(k,1); m = terms(k,2);
    R = radial_poly(n, abs(m), rho);
    if m > 0
        B(:,k) = R .* cos(m * theta);
    elseif m < 0
        B(:,k) = R .* sin(abs(m) * theta);
    else
        B(:,k) = R;
    end
end
end

function terms = zernike_term_list(nMax)
terms = [];
for n = 0:nMax
    for m = -n:2:n
        terms = [terms; n, m]; %#ok<AGROW>
    end
end
end

function R = radial_poly(n, m, rho)
R = zeros(size(rho));
for s = 0:((n-m)/2)
    c = (-1)^s * factorial(n-s) / (factorial(s) * factorial((n+m)/2-s) * factorial((n-m)/2-s));
    R = R + c * rho.^(n - 2*s);
end
end
