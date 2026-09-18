function [masterDisplacement,info] = MiuraPerturbCorrect( ...
    targetMasterDisplacement,constraintDirections,referenceComponents, ...
    gij,options)
%MIURAPERTURBCORRECT Project a perturbed state onto fixed quadratic constraints.
%
% The solved proximal problem is
%   min 0.5*||c(deltaXm)||^2
%       + 0.5*lambda*||deltaXm-targetDeltaXm||^2,
% using the analytic Jacobian, exact residual curvature, and adaptive
% Levenberg-Marquardt damping. No chirality or intersection constraints are
% evaluated in this first kinematic version.

if isempty(options)
    options = struct();
end
regularization = optionValue(options,'Regularization',1e-6);
constraintTolerance = optionValue(options,'ConstraintTolerance',1e-6);
stepTolerance = optionValue(options,'StepTolerance',1e-11);
gradientTolerance = optionValue(options,'GradientTolerance',1e-10);
maxIterations = optionValue(options,'MaxIterations',80);
maxInnerIterations = optionValue(options,'MaxInnerIterations',25);
verbose = optionValue(options,'Verbose',false);

if regularization <= 0 || constraintTolerance <= 0 || ...
        stepTolerance <= 0 || gradientTolerance <= 0 || ...
        maxIterations < 1 || maxInnerIterations < 1
    error('MiuraPerturbCorrect:InvalidOptions', ...
        'Solver tolerances, regularization, and iteration limits must be positive.');
end

targetMasterDisplacement = targetMasterDisplacement(:);
masterDisplacement = targetMasterDisplacement;
% Validate once and cache fixed factors outside the iterative hot path.
[residual,jacobian,residualCurvature] = AffineMetricResidualJacobian( ...
    masterDisplacement,constraintDirections,referenceComponents,gij);
directionU = constraintDirections(:,:,1);
directionV = constraintDirections(:,:,2);
referenceU = referenceComponents(:,:,1);
referenceV = referenceComponents(:,:,2);
gij = gij(:);
nMasterDof = numel(masterDisplacement);
diagonalIndices = 1:nMasterDof+1:nMasterDof^2;
mu = [];
exitflag = 0;
iterations = 0;

for iteration = 1:maxIterations
    iterations = iteration;

    % Use exact first and second derivatives of the quadratic constraints.
    if iteration > 1
        [residual,jacobian,residualCurvature] = EvaluateLowRankConstraints( ...
            masterDisplacement,directionU,directionV,referenceU,referenceV,gij);
    end
    targetDifference = masterDisplacement-targetMasterDisplacement;
    objective = 0.5*(residual.'*residual) ...
        +0.5*regularization*(targetDifference.'*targetDifference);
    gradient = jacobian.'*residual+regularization*targetDifference;
    hessian = jacobian.'*jacobian+residualCurvature;
    hessian(diagonalIndices) = hessian(diagonalIndices)+regularization;
    hessian = 0.5*(hessian+hessian.');

    if any(~isfinite([objective;gradient;hessian(:)]))
        exitflag = -2;
        break
    end

    if norm(gradient,inf) <= gradientTolerance
        exitflag = 1;
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
        if norm(step) <= stepTolerance*displacementScale
            exitflag = 3;
            break
        end

        trialMasterDisplacement = masterDisplacement+step;
        residualTrial = EvaluateLowRankConstraints( ...
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

% Report objective stationarity and feasibility as separate diagnostics.
[finalResidual,finalJacobian] = EvaluateLowRankConstraints( ...
    masterDisplacement,directionU,directionV,referenceU,referenceV,gij);
finalTargetDifference = masterDisplacement-targetMasterDisplacement;
finalGradient = finalJacobian.'*finalResidual ...
    +regularization*finalTargetDifference;
maxResidual = norm(finalResidual,inf);
gradientInfinityNorm = norm(finalGradient,inf);
if exitflag == 0 && gradientInfinityNorm <= gradientTolerance
    exitflag = 1;
elseif exitflag == 3
    exitflag = 2;
end

info = struct();
info.exitflag = exitflag;
info.iterations = iterations;
info.finalDamping = mu;
info.maxResidual = maxResidual;
info.gradientInfinityNorm = gradientInfinityNorm;
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
