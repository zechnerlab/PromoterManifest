function conditions = ConvertConcNaming()
%% DM case

time = [10, 20, 30, 40, 50];
conc = [100, 275, 690, 3000];

conditions = {};
counter = 1;
for i=1:length(time)
    
    for k=1:length(conc)
        condition.PulseParameters = [1, time(i)];
        condition.Concentration = conc(k);
        if (condition.Concentration < 1000)
            condition.Name = ['DM_' num2str(time(i)) 'min_' num2str(conc(k)) 'nM'];
        else
            condition.Name = ['DM_' num2str(time(i)) 'min_' num2str(floor(conc(k)/1000)) 'uM'];
        end
        conditions{counter} = condition;
        counter = counter + 1;
    end
    
end

%% FM case

numPulses = [2, 3, 4, 5, 6, 8];
conc = [690];
time = 5;

for i=1:length(numPulses)
    
    for k=1:length(conc)
        condition.PulseParameters = [numPulses(i), time, time];
        condition.Concentration = conc(k);
        
        if (condition.Concentration < 1000)
            condition.Name = ['FM_' num2str(numPulses(i)) '_' num2str(time) 'min_' num2str(conc(k)) 'nM'];
        else
            condition.Name = ['FM_' num2str(numPulses(i)) '_' num2str(floor(conc(k)/1000)) 'uM'];
        end
        conditions{counter} = condition;
        counter = counter + 1;
    end
    
end

end