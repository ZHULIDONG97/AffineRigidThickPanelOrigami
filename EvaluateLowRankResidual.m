function [residual,projectedU,projectedV] = EvaluateLowRankResidual( ...
    beta,directionU,directionV,referenceU,referenceV,gij)
% Evaluate prevalidated constraints, retaining projections for an accepted trial.
betaByCoordinate = reshape(beta,size(directionU,2),3);
projectedU = directionU*betaByCoordinate;
projectedV = directionV*betaByCoordinate;
residual = sum(projectedU.*projectedV ...
    +referenceU.*projectedV+referenceV.*projectedU,2)-gij;
end
