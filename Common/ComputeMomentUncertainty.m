function [Moments, Uncertainties, NumCells]=ComputeMomentUncertainty(data)

    timePoints = size(data, 2);

    Moments = zeros(2, timePoints);
    Uncertainties = zeros(2, timePoints);

    for k=1:timePoints
        N = length(data{k});

        firstMoment = mean(data{k});
        secondMoment = 1/N*sum(power(data{k},2));


        fourthCentralMoment = 1/N * sum(power(data{k} - repmat(firstMoment, N, 1), 4));

        Moments(1, k) = firstMoment;
        Moments(2, k) = secondMoment - firstMoment.^2;

        Uncertainties(1, k) = 1/N * Moments(2, k);
        Uncertainties(2, k) = 1/N * (fourthCentralMoment - ...
            (N - 3)/(N - 1).*Moments(2, k).^2);

            
        NumCells(k) = N;
    end

end