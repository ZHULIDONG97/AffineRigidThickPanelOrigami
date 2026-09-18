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
    r = size(directionU,2);
    % Reuse the same reduced curvature block for the three coordinates.
    coordinateCurvature = directionU.'*(residual.*directionV) ...
        +directionV.'*(residual.*directionU);
    residualCurvature = zeros(3*r,3*r);
    for coordinate = 1:3
        indices = (coordinate-1)*r+(1:r);
        residualCurvature(indices,indices) = coordinateCurvature;
    end
end
end
