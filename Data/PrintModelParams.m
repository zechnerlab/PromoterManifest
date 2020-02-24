function s = PrintModelParams(model)

    s = sprintf('&$%.2e$ &$%.2e$ &$%.2e$ &$%.2e$ &$%.2e$ &$%.2e$ &$%.2e$ &$%.2e$ &$%.2e$ &$%.2e$ &$%.2e$',...
        model.Q(2, 1), model.Q(1, 2), model.Q(3, 2), model.Q(2, 3), model.Z(2), model.Z(3), ...
        model.c(2), model.c(4), model.c(3)*model.H(1), model.c(3)^2*model.H(2), model.MeasurementSigma);

end