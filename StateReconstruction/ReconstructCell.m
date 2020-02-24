function particles = ReconstructCell(model, config, particles, cellStruct, filterLength)
    
    
   M = length(particles);
   numMeasurements = length(cellStruct.Measurement);
   numRecursions = ceil(numMeasurements / (filterLength-1))-1;
   measurementTimes = cellStruct.MeasurementTime;

    
    LAll = ones(M, 1);
    
    for l=1:numRecursions
        
        measureIdx = ((l-1)*(filterLength-1) + 1):(min(l*(filterLength-1)+1, numMeasurements));
        nextT = measurementTimes(measureIdx(end));
        
        
        LAll(or(isnan(LAll), isinf(LAll))) = -inf;
        
        if (sum(abs(imag(LAll))>0)>0)
            LAll(abs(imag(LAll))>0) = -inf;
        end
        
        W = exp(LAll - max(LAll));
        W = W/sum(W);
        if (sum(isnan(W))>0 || sum(W<0)>0)
            particles = {};
           return; 
        end
        idx = randsample(1:M, M, true, W);
        particles = particles(idx);
        
        [particles] = SimulatePromoter(model, particles, nextT);
        
%         figure(98);
%         for k=1:length(particles)
%            stairs(particles{k}.t/60, particles{k}.Z); hold on;
%         end
%         
%         hold off;
%         
        
        for k=1:M
            
           
            
            M0 = particles{k}.PosteriorStats(:, end);
            
            [Mu, Sigma, PosteriorTime, Likelihood, PosteriorStats] = FilterInputMoments(particles{k}.Model, config, particles{k}.t, particles{k}.Z, cellStruct, M0, measureIdx);
            
            
            particles{k}.PosteriorStats = [particles{k}.PosteriorStats, PosteriorStats(:, 2:end)];
            particles{k}.PosteriorTime = [particles{k}.PosteriorTime, PosteriorTime(2:end)];
            particles{k}.L = sum(Likelihood);
            LAll(k) = particles{k}.L;
            
%             subplot(1,2,1);
%             stairs(particles{k}.t, particles{k}.Z);
%             
%             subplot(1,2,2);
%             plot(particles{k}.PosteriorTime/60, particles{k}.PosteriorStats(2, :), cellStruct.MeasurementTime/60, cellStruct.Measurement, 'o');
%             drawnow;
            
            %fprintf('%d\n', k);
        end
        
        W = exp(LAll - max(LAll));
        W = W/sum(W);
        
        for k=1:M
            particles{k}.W = W(k);
        end
    end
    
    
    
    