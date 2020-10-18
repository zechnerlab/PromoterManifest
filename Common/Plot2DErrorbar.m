function p = Plot2DErrorbar(data1, data2, data1sem, data2sem, color, semFactor)

for j=1:length(data1)
    p(1) = plot(data1(j) + semFactor*[-data1sem(j), data1sem(j)], [data2(j), data2(j)]);
    set(p, 'Color', color, 'LineWidth', 1.5);
    p(2) = plot([data1(j), data1(j)], data2(j) + semFactor*[-data2sem(j), data2sem(j)]);
    set(p, 'Color', color, 'LineWidth', 1.5);

end