function X0 = SetInitialInputMoments(infos, X0, means, vars)

secondMoments = means.^2 + vars;

k = [];

    % set means and second order moments of single species
    for l=1:length(means)
        X0 = SetMoment(infos.MomentSystem, X0, k, l, means(l));
        X0 = SetMoment(infos.MomentSystem, X0, k, [l l], secondMoments(l));
    end
    
    % set cross moments
    
    X0 = SetMoment(infos.MomentSystem, X0, k, [1, 2], means(1)*means(2));
    X0 = SetMoment(infos.MomentSystem, X0, k, [1, 3], means(1)*means(3));
%    X0 = SetMoment(infos.MomentSystem, X0, k, [1, 4], means(1)*means(4)*X0(k));
    X0 = SetMoment(infos.MomentSystem, X0, k, [2, 3], means(2)*means(3));
 %   X0 = SetMoment(infos.MomentSystem, X0, k, [2, 4], means(2)*means(4)*X0(k));
 %   X0 = SetMoment(infos.MomentSystem, X0, k, [3, 3], means(3)*means(3)*X0(k));
  %  X0 = SetMoment(infos.MomentSystem, X0, k, [3, 4], means(3)*means(4)*X0(k));
    
    X0 = SetMoment(infos.MomentSystem, X0, k, [1, 1, 3], secondMoments(1)*means(3));
    X0 = SetMoment(infos.MomentSystem, X0, k, [1, 2, 3], means(1)*means(2)*means(3));
    X0 = SetMoment(infos.MomentSystem, X0, k, [1, 3, 3], means(1)*secondMoments(3));
    
    X0 = SetMoment(infos.MomentSystem, X0, k, [1, 1, 3, 3], secondMoments(1)*secondMoments(3));

