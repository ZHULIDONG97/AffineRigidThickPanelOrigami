function [masterDisplacement,info] = MiuraPerturbCorrect( ...
    targetMasterDisplacement,constraintDirections,referenceComponents, ...
    gij,options)
%MIURAPERTURBCORRECT Project a perturbed state onto fixed quadratic constraints.
%
% The solved proximal problem is
%   min 0.5*||c(deltaXm)||^2
%       + 0.5*lambda*||deltaXm-targetDeltaXm||^2,
% using the analytic Jacobian, residual curvature, and adaptive
% Levenberg-Marquardt damping. No chirality or intersection constraints are
% evaluated in this first kinematic version.
% Constant curvature coefficients below a relative 1e-12 cutoff are pruned.
% Stopping matches RigidOrigamiSimulator: scaled RMS residual <= 1e-8,
% scaled gradient 2-norm <= 1e-6, or relative step <= 1e-6 (exitflags 1/2/3).
% ResidualTolerance, GradientTolerance and StepTolerance override these
% defaults. ConstraintTolerance independently tests the raw maximum residual.

if isempty(options)
    options = struct();
end
regularization = optionValue(options,'Regularization',1e-6);
constraintTolerance = optionValue(options,'ConstraintTolerance',1e-6);
residualTolerance = optionValue(options,'ResidualTolerance',1e-8);
stepTolerance = optionValue(options,'StepTolerance',1e-6);
gradientTolerance = optionValue(options,'GradientTolerance',1e-6);
maxIterations = optionValue(options,'MaxIterations',80);
maxInnerIterations = optionValue(options,'MaxInnerIterations',25);
verbose = optionValue(options,'Verbose',false);

if regularization <= 0 || constraintTolerance <= 0 || residualTolerance <= 0 || ...
        stepTolerance <= 0 || gradientTolerance <= 0 || ...
        maxIterations < 1 || maxInnerIterations < 1
    error('MiuraPerturbCorrect:InvalidOptions', ...
        'Solver tolerances, regularization, and iteration limits must be positive.');
end

targetMasterDisplacement = targetMasterDisplacement(:);
masterDisplacement = targetMasterDisplacement;
% Validate once and cache fixed factors outside the iterative hot path.
[residual,jacobian] = AffineMetricResidualJacobian( ...
    masterDisplacement,constraintDirections,referenceComponents,gij);
directionU = constraintDirections(:,:,1);
directionV = constraintDirections(:,:,2);
referenceU = referenceComponents(:,:,1);
referenceV = referenceComponents(:,:,2);
gij = gij(:);
% Match the geometric scaling and stopping measures in RigidOrigamiSimulator.
nConstraints = numel(gij);
referenceMagnitude = sqrt(sum(referenceU.^2,2).*sum(referenceV.^2,2));
constraintScale = max([1;referenceMagnitude;abs(gij)]);
nMasterDof = numel(masterDisplacement);
diagonalIndices = 1:nMasterDof+1:nMasterDof^2;
% Cache G(:,i)=Hi(:) across calls; only the constraint directions determine Hi.
persistent cachedDirections curvatureFactors
if ~isequal(cachedDirections,constraintDirections)
    nFreeMasters = size(directionU,2);
    [curvatureRow,curvatureColumn] = ndgrid(1:nFreeMasters);
    coordinateCurvature = (directionU(:,curvatureRow(:)).*directionV(:,curvatureColumn(:)) ...
        +directionV(:,curvatureRow(:)).*directionU(:,curvatureColumn(:))).';
    % Remove tiny coefficients relative to the model's largest Hi entry.
    curvatureTolerance = 1e-12*max(abs(coordinateCurvature(:)));
    if isfinite(curvatureTolerance)
        coordinateCurvature(abs(coordinateCurvature)<=curvatureTolerance) = 0;
    end
    [blockIndex,constraintIndex,curvatureValue] = find(coordinateCurvature);
    hessianIndex = curvatureRow(blockIndex) ...
        +(curvatureColumn(blockIndex)-1)*nMasterDof;
    hessianIndex = hessianIndex(:);
    blockOffset = nFreeMasters*(nMasterDof+1);
    % Hi repeats its coordinate block three times, with zero cross blocks.
    curvatureFactors = sparse([hessianIndex;hessianIndex+blockOffset; ...
        hessianIndex+2*blockOffset],repmat(constraintIndex,3,1), ...
        repmat(curvatureValue,3,1),nMasterDof^2,nConstraints);
    cachedDirections = constraintDirections;
end
mu = [];
exitflag = 0;
iterations = 0;
relativeStep = inf;

for iteration = 1:maxIterations
    iterations = iteration;

    % Residual and Jacobian belong to the current accepted state.
    targetDifference = masterDisplacement-targetMasterDisplacement;
    objective = 0.5*(residual.'*residual) ...
        +0.5*regularization*(targetDifference.'*targetDifference);
    gradient = jacobian.'*residual+regularization*targetDifference;
    if any(~isfinite([objective;gradient]))
        exitflag = -2;
        break
    end

    constraintMeasure = norm(residual)/(sqrt(nConstraints)*constraintScale);
    gradientMeasure = norm(gradient)/max(1,sqrt(2*objective));
    if constraintMeasure <= residualTolerance
        exitflag = 1;
        break
    end
    if gradientMeasure <= gradientTolerance
        exitflag = 2;
        break
    end

    % Accumulate all ci*Hi in one sparse multiply; retain the dense small solve.
    hessian = jacobian.'*jacobian ...
        +reshape(full(curvatureFactors*residual),nMasterDof,nMasterDof);
    hessian(diagonalIndices) = hessian(diagonalIndices)+regularization;
    hessian = 0.5*(hessian+hessian.');
    if any(~isfinite(hessian(:)))
        exitflag = -2;
        break
    end

    % Initialize damping above any negative curvature in the exact Hessian.
    if isempty(mu)
        hessianScale = max(1,norm(hessian,inf));
        minimumEigenvalue = min(eig(hessian));
        mu = max(1e-8*hessianScale, ...
            -minimumEigenvalue+sqrt(eps)*hessianScale);
        mu = max(mu,eps*hessianScale);
        muMaximum = 1e12*hessianScale;
    end

    accepted = false;
    displacementScale = max(1,norm(masterDisplacement));
    for innerIteration = 1:maxInnerIterations
        % Cholesky both detects and solves a positive-definite damped system.
        dampedHessian = hessian;
        dampedHessian(diagonalIndices) = dampedHessian(diagonalIndices)+mu;
        [R,notPositiveDefinite] = chol(dampedHessian);
        if notPositiveDefinite ~= 0
            mu = min(2*mu,muMaximum);
            continue
        end
        step = -R\(R.'\gradient);
        relativeStep = norm(step)/displacementScale;
        if relativeStep <= stepTolerance
            exitflag = 3;
            break
        end

        trialMasterDisplacement = masterDisplacement+step;
        [residualTrial,projectedUTrial,projectedVTrial] = EvaluateLowRankResidual( ...
            trialMasterDisplacement,directionU,directionV,referenceU,referenceV,gij);
        trialTargetDifference = trialMasterDisplacement ...
            -targetMasterDisplacement;
        objectiveTrial = 0.5*(residualTrial.'*residualTrial) ...
            +0.5*regularization ...
            *(trialTargetDifference.'*trialTargetDifference);
        actualReduction = objective-objectiveTrial;
        predictedReduction = -(gradient.'*step ...
            +0.5*step.'*hessian*step);
        if predictedReduction > 0
            reductionRatio = actualReduction/predictedReduction;
        else
            reductionRatio = -inf;
        end

        if reductionRatio > 1e-4 && actualReduction > 0
            masterDisplacement = trialMasterDisplacement;
            % Reuse the accepted residual/projections; rejected trials never overwrite them.
            residual = residualTrial;
            jacobian = EvaluateLowRankJacobian(projectedUTrial,projectedVTrial, ...
                directionU,directionV,referenceU,referenceV);
            mu = max(mu*max(1/3,1-(2*reductionRatio-1)^3), ...
                eps*max(1,norm(hessian,inf)));
            accepted = true;
            break
        end
        mu = min(2*mu,muMaximum);
    end

    if exitflag == 3 || ~accepted
        if exitflag == 0
            exitflag = -1;
        end
        break
    end
end

% Reuse current-state derivatives, including the last accepted iteration-limit step.
finalTargetDifference = masterDisplacement-targetMasterDisplacement;
finalGradient = jacobian.'*residual ...
    +regularization*finalTargetDifference;
maxResidual = norm(residual,inf);
gradientInfinityNorm = norm(finalGradient,inf);
finalObjective = 0.5*(residual.'*residual) ...
    +0.5*regularization*(finalTargetDifference.'*finalTargetDifference);
if isempty(mu)
    mu = 0;
end

info = struct();
info.exitflag = exitflag;
info.iterations = iterations;
info.finalDamping = mu;
info.maxResidual = maxResidual;
info.gradientInfinityNorm = gradientInfinityNorm;
info.constraintMeasure = norm(residual)/(sqrt(nConstraints)*constraintScale);
info.gradientMeasure = norm(finalGradient)/max(1,sqrt(2*finalObjective));
info.relativeStep = relativeStep;
info.objective = finalObjective;
% Solver termination and the caller's absolute feasibility test are independent.
info.constraintSatisfied = maxResidual <= constraintTolerance;
if verbose
    fprintf(['Perturb-correct: exit %d, iterations %d, ' ...
        'max |c| %.3e, max |grad| %.3e, damping %.3e.\n'], ...
        exitflag,iterations,maxResidual,info.gradientInfinityNorm,mu);
end
end

function value = optionValue(options,name,defaultValue)
if isfield(options,name)
    value = options.(name);
else
    value = defaultValue;
end
end
