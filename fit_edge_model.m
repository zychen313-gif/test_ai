function model = fit_edge_model(trainData, signalType, strategy, cfg)
[X, Y, R, T, Z] = gather_train(trainData, signalType);

switch strategy
    case 'mean_previous'
        model.X = X; model.Y = Y; model.Z = Z;

    case 'rbf_gaussian'
        n = numel(Z);
        m = min(cfg.rbfCenters, n);
        idx = round(linspace(1, n, m));
        C = [X(idx), Y(idx)];
        Phi = rbf_kernel([X,Y], C, cfg.rbfGamma);
        w = (Phi' * Phi + cfg.ridgeLambda * eye(size(Phi,2))) \ (Phi' * Z);
        model.C = C; model.w = w; model.gamma = cfg.rbfGamma;

    case 'bspline_tensor'
        Br = bspline_design(R, cfg.rMin, cfg.rMax, cfg.bsplineOrder, cfg.bsplineRadialKnots);
        Bt = bspline_design(T, 0, 360, cfg.bsplineOrder, cfg.bsplineThetaKnots);
        D = kron_columns(Br, Bt);
        w = (D' * D + cfg.ridgeLambda * eye(size(D,2))) \ (D' * Z);
        model.w = w; model.meta = [size(Br,2), size(Bt,2)];
        model.order = cfg.bsplineOrder;
        model.nr = cfg.bsplineRadialKnots;
        model.nt = cfg.bsplineThetaKnots;
        model.rMin = cfg.rMin;
        model.rMax = cfg.rMax;

    case 'zernike'
        B = zernike_design(R, T, cfg.rMax, cfg.zernikeOrder);
        w = (B' * B + cfg.ridgeLambda * eye(size(B,2))) \ (B' * Z);
        model.w = w; model.zOrder = cfg.zernikeOrder;
        model.rMax = cfg.rMax;

    case 'fourier_radial'
        B = fourier_radial_design(R, T);
        w = (B' * B + cfg.ridgeLambda * eye(size(B,2))) \ (B' * Z);
        model.w = w;

    otherwise
        error('Unsupported strategy: %s', strategy);
end
end

function [X, Y, R, T, Z] = gather_train(trainData, signalType)
X=[];Y=[];R=[];T=[];Z=[];
for i=1:numel(trainData)
    d = trainData{i};
    m = d.edgeMask;
    X = [X; d.x(m)]; %#ok<AGROW>
    Y = [Y; d.y(m)]; %#ok<AGROW>
    R = [R; d.r(m)]; %#ok<AGROW>
    T = [T; d.theta(m)]; %#ok<AGROW>
    switch signalType
        case 'X', Zi = d.zX(m);
        case 'Y', Zi = d.zY(m);
    end
    Z = [Z; Zi]; %#ok<AGROW>
end
end

function Phi = rbf_kernel(P, C, gamma)
D2 = sqdist(P, C);
Phi = exp(-gamma * D2);
end


function D2 = sqdist(A, B)
AA = sum(A.^2,2);
BB = sum(B.^2,2)';
D2 = max(AA + BB - 2*(A*B'), 0);
end
