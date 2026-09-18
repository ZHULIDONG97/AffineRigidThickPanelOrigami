function [residual,jacobian,residualCurvature] = ...
    AffineMetricResidualJacobian(beta,constraintDirections, ...
    referenceComponents,gij)
%AFFINEMETRICRESIDUALJACOBIAN Evaluate selected A'*A-I constraints exactly.

% Validate the cached low-rank factors before evaluating the polynomial.
beta = beta(:);
nc = size(constraintDirections,1);
r = size(constraintDirections,2);
if size(constraintDirections,3) ~= 2 || ...
        ~isequal(size(referenceComponents),[nc,3,2]) || ...
        numel(gij) ~= nc || numel(beta) ~= 3*r
    error('AffineMetricResidualJacobian:InvalidConstraintData', ...
        'Constraint factors and the coordinate-block beta vector are inconsistent.');
end
if any(~isfinite(beta)) || any(~isfinite(constraintDirections(:))) || ...
        any(~isfinite(referenceComponents(:))) || any(~isfinite(gij))
    error('AffineMetricResidualJacobian:NonfiniteInput', ...
        'All metric-constraint inputs must be finite.');
end

directionU = constraintDirections(:,:,1);
directionV = constraintDirections(:,:,2);
referenceU = referenceComponents(:,:,1);
referenceV = referenceComponents(:,:,2);
betaByCoordinate = reshape(beta,r,3);
projectedU = directionU*betaByCoordinate;
projectedV = directionV*betaByCoordinate;

% Evaluate (a_u+du)'*(a_v+dv)-delta without assembling a large Hessian.
residual = sum(projectedU.*projectedV ...
    +referenceU.*projectedV+referenceV.*projectedU,2)-gij(:);

if nargout > 1
    % The Jacobian is linear in beta for every fixed quadratic constraint.
    coefficientU = projectedV+referenceV;
    coefficientV = projectedU+referenceU;
    jacobian = [ ...
        coefficientU(:,1).*directionU+coefficientV(:,1).*directionV, ...
        coefficientU(:,2).*directionU+coefficientV(:,2).*directionV, ...
        coefficientU(:,3).*directionU+coefficientV(:,3).*directionV];
end

if nargout > 2
    % Sum residual_k*Hessian(residual_k) for the exact Newton Hessian.
    coordinateCurvature = directionU.'*(residual.*directionV) ...
        +directionV.'*(residual.*directionU);
    residualCurvature = zeros(3*r,3*r);
    for coordinate = 1:3
        indices = (coordinate-1)*r+(1:r);
        residualCurvature(indices,indices) = coordinateCurvature;
    end
end
end
