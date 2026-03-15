function cfg = default_config()
cfg.rMin = 100;
cfg.rMax = 150;
cfg.notchCenterDeg = 270;
cfg.notchHalfWidthDeg = 10;
cfg.outlierSigma = 3;
cfg.gridN = 241;
cfg.minTrainWafers = 1;
cfg.ridgeLambda = 1e-4;
cfg.rbfGamma = 2.5e-3;
cfg.rbfCenters = 400;
cfg.bsplineOrder = 4;
cfg.bsplineRadialKnots = 10;
cfg.bsplineThetaKnots = 16;
cfg.zernikeOrder = 6;
cfg.strategies = {'mean_previous','rbf_gaussian','bspline_tensor','zernike','fourier_radial'};
end
