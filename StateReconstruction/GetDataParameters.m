function [promoters, concentrations, numPulses, durations, intervals] = GetDataParameters(dataB)

    promoters = unique(dataB.Promoter)';
    concentrations = unique(dataB.Condition(:, 1));
    numPulses = unique(dataB.Condition(:, 2));
    durations = unique(dataB.Condition(:, 3));
    intervals = unique(dataB.Condition(:, 4));
    intervals = setdiff(intervals, -1);
end