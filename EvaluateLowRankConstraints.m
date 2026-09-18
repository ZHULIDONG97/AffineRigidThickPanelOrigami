function [residual,jacobian,residualCurvature] = EvaluateLowRankConstraints( ...
    beta,directionU,directionV,referenceU,referenceV,gij)
% Evaluate prevalidated quadratic factors, as in RigidOrigamiSimulator's LM.

r = size(directionU,2);
betaByCoordinate = reshape(beta,r,3);

% Preserve the original product order while reusing the fixed factor rows.
projectedU = directionU*betaByCoordinate;
projectedV = directionV*betaByCoordinate;
residual = sum(projectedU.*projectedV ...
    +referenceU.*projectedV+referenceV.*projectedU,2)-gij;

if nargout > 1
    % Each quadratic constraint has a Jacobian linear in the displacement.
    coefficientU = projectedV+referenceV;
    coefficientV = projectedU+referenceU;
    jacobian = [ ...
        coefficientU(:,1).*directionU+coefficientV(:,1).*directionV, ...
        coefficientU(:,2).*directionU+coefficientV(:,2).*directionV, ...
        coefficientU(:,3).*directionU+coefficientV(:,3).*directionV];
end

if nargout > 2
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
