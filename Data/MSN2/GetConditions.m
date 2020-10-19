function conditionsOut = GetConditions(conditions, pulseLength, conc)
%% DM case

counter = 1;
for k=1:length(conditions)
   
    if (sum(conditions{k}.PulseParameters(2) == pulseLength)>0)
       if  (sum(conditions{k}.Concentration == conc)>0)
          conditionsOut{counter} = conditions{k};
          counter = counter + 1;
       end
    end
    
end