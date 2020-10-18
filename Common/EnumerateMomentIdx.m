function moments = EnumerateMomentIdx(numSpecies, orderVec)

counter = 1;

for k=1:length(orderVec)
    order = orderVec(k);
    PermMat = npermutek(1:numSpecies, order);
    PermMat = sort(PermMat, 2);
    momentIdx = unique(PermMat, 'rows');
    
    for u=1:size(momentIdx, 1)
        moments{counter}.Species = momentIdx(u, :);
        counter = counter + 1;
    end
end

end