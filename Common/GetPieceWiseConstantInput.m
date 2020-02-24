function [u] = GetPieceWiseConstantInput(t, params)

    inputLevels = params.inputLevels;
    inputTimes = params.inputTimes;
    
    tIdx = find(t<=inputTimes);
    
    %tIdx = sum(t>inputTimes)+1;
    
    if (isempty(tIdx))
       u = inputLevels(1); 
    else
       u = inputLevels(tIdx(1)); 
    end
    
    
%     for l=1:length(inputTimes)
%        if t <= inputTimes(l)
%            u = inputLevels(l);
%            return;
%        end
%     end
%   
%     u = inputLevels(1);
    
%     val = inputLevels(inputTimes >= t);
%     
%     if (isempty(val))
%        u = inputLevels(1); 
%     else
%        u = val(1);
%     end
%     
%     u = abs(u);
end