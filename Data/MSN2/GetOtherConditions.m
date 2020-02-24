function conditionsOut = GetOtherConditions(conditions, excludeConditions)
%% DM case

counter = 1;
for k=1:length(conditions)
    include = 1;
    for u=1:length(excludeConditions)
        if (strcmp(conditions{k}.Name, excludeConditions{u}.Name)==1)
            include = 0;
            break;
        end
    end
    
    if (include)
        conditionsOut{counter} = conditions{k};
        counter = counter + 1;
    end
end