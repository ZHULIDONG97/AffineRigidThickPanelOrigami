function [residual,jacobian,residualCurvature] = EvaluateLowRankConstraints( ...
    beta,directionU,directionV,referenceU,referenceV,gij)
% Evaluate prevalidated quadratic factors, as in RigidOrigamiSimulator's LM.

[residual,projectedU,projectedV] = EvaluateLowRankResidual( ...
    beta,directionU,directionV,referenceU,referenceV,gij);

if nargout > 1
    jacobian = EvaluateLowRankJacobian( ...
        projectedU,projectedV,directionU,directionV,referenceU,referenceV);
end

if nargout > 2
    % Reuse the same reduced curvature block for the three coordinates.
    coordinateCurvature = directionU.'*(residual.*directionV) ...
        +directionV.'*(residual.*directionU);
    residualCurvature = blkdiag(coordinateCurvature,coordinateCurvature,coordinateCurvature);
end
end
