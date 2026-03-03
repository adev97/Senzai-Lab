function params = fit2DGaussian(rfMap)
    % rfMap = 2D receptive field (matrix)
    [X,Y] = meshgrid(1:size(rfMap,2), 1:size(rfMap,1));
    xdata = [X(:), Y(:)];
    zdata = rfMap(:);

    % Initial guess: amplitude, x0, y0, sigma_x, sigma_y, offset
    A0 = max(zdata);
    x0 = size(rfMap,2)/2;
    y0 = size(rfMap,1)/2;
    sigma_x0 = size(rfMap,2)/4;
    sigma_y0 = size(rfMap,1)/4;
    offset = min(zdata);

    % Gaussian function
    gauss2D = @(p,x) p(1)*exp(-(((x(:,1)-p(2)).^2)/(2*p(4)^2) + ((x(:,2)-p(3)).^2)/(2*p(5)^2))) + p(6);

    opts = optimoptions('lsqcurvefit','Display','off');
    lb = [0, 1, 1, 0.1, 0.1, 0];
    ub = [Inf, size(rfMap,2), size(rfMap,1), Inf, Inf, max(zdata)];
    pfit = lsqcurvefit(gauss2D, [A0,x0,y0,sigma_x0,sigma_y0,offset], xdata, zdata, lb, ub, opts);

    params.A       = pfit(1);
    params.x0      = pfit(2);
    params.y0      = pfit(3);
    params.sigma_x = pfit(4);
    params.sigma_y = pfit(5);
    params.offset  = pfit(6);
end