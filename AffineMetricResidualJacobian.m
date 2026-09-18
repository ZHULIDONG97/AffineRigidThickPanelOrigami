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

% Keep validation at the public boundary; LM caches these fixed factors once.
directionU = constraintDirections(:,:,1);
directionV = constraintDirections(:,:,2);
referenceU = referenceComponents(:,:,1);
referenceV = referenceComponents(:,:,2);
if nargout > 2
    [residual,jacobian,residualCurvature] = EvaluateLowRankConstraints( ...
        beta,directionU,directionV,referenceU,referenceV,gij(:));
elseif nargout > 1
    [residual,jacobian] = EvaluateLowRankConstraints( ...
        beta,directionU,directionV,referenceU,referenceV,gij(:));
else
    residual = EvaluateLowRankConstraints( ...
        beta,directionU,directionV,referenceU,referenceV,gij(:));
end
end
