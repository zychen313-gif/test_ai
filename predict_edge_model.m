function zPred = predict_edge_model(model, x, y, r, t, strategy)
switch strategy
    case 'mean_previous'
        F = scatteredInterpolant(model.X, model.Y, model.Z, 'natural', 'nearest');
        zPred = F(x, y);

    case 'rbf_gaussian'
        D2 = sqdist([x,y], model.C);
        Phi = exp(-model.gamma * D2);
        zPred = Phi * model.w;

    case 'bspline_tensor'
        Br = bspline_design(r, model.rMin, model.rMax, model.order, model.nr);
        Bt = bspline_design(t, 0, 360, model.order, model.nt);
        D = kron_columns(Br, Bt);
        zPred = D * model.w;

    case 'zernike'
        B = zernike_design(r, t, model.rMax, model.zOrder);
        zPred = B * model.w;

    case 'fourier_radial'
        B = fourier_radial_design(r, t);
        zPred = B * model.w;

    otherwise
        error('Unsupported strategy: %s', strategy);
end
end


function D2 = sqdist(A, B)
AA = sum(A.^2,2);
BB = sum(B.^2,2)';
D2 = max(AA + BB - 2*(A*B'), 0);
end
