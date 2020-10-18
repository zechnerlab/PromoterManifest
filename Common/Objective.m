function [L, postStats, postTimes, model] = Objective(commonParams, model, config, targetData)

    if sum(commonParams<0)>0
       L = -inf;
       postStats = 0;
       postTimes = 0;
       model = 0;
       return;
    end

    model = SetModelParameters(model, config, commonParams);
    
    if sum(model.P0<0)>0
       L = -inf;
       fprintf('Invalid initial condition\n');
       postStats = 0;
       postTimes = 0;
       model = 0;
       return;
    end
        
    cells = targetData.Conditions{1}.cells;
    parfor k=1:length(cells)
        warning off;
        options = model;
        [PosteriorStats, PosteriorTime, LTrace] = RunFilter(cells{k}, options);

    
        postStats{k} = PosteriorStats;
        postTimes{k} = PosteriorTime;
        
        LVec(k) = sum(sum(real(LTrace)));
        
        %fprintf('Evaluated trace %d...\n', k);
        
    end
    L = sum(LVec) + EvaluatePrior(model, config);
    
    if (~isreal(L))
       L = -inf; 
    end
    
    if (isinf(L))
       L = -inf; 
    end
    
%     figure(77);
%     counter = 0;
%     K = 2;
%     numStates = size(model.Q, 1);
%     for k=1:K%length(cells)
%         postStat = postStats{k};
%         postTime = postTimes{k};
%         
%         counter = counter + 1;
%         subplot(K,2,counter);
%         [v, MAPState] = max(postStat(1:numStates, :));
%         stairs(postTime, postStat(1, :)); hold on;
%         stairs(postTime, postStat(2, :));
%         stairs(postTime, postStat(3, :)); hold off;
%         %stairs(postTime, postStat(4, :)); hold off;
%         
%         counter = counter + 1;
%         subplot(K,2,counter);
%         plot(cells{k}.MeasurementTime, cells{k}.Measurement, 'o');
%         
%     end
%     drawnow;

end