function [strOut, momentName]= SubstituteMoment(names, strIn, maxOrder)

    N = length(names);
    idx = 1:N;
    
    strOut = strIn;
    counter = 1;
    
    for k=maxOrder:-1:1
        [PermMat, S] = npermutek(idx, k);
        
        for i=1:size(PermMat, 1)
           
            searchStr = char(names(PermMat(i, 1), :));
            replaceStr = char(names(PermMat(i, 1), :));
            for j=2:k
               searchStr = sprintf('%s*%s', searchStr, char(names(PermMat(i, j), :)));
               replaceStr = sprintf('%s%d', replaceStr, PermMat(i, j));
            end
            
            strOut = subs(strOut, simplify(sym(searchStr)), replaceStr);
            momentName{counter} = replaceStr;
            counter = counter + 1;
        end
    end

end