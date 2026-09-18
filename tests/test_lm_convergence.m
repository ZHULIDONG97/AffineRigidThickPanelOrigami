function tests = test_lm_convergence
% Verify the reference LM stopping measures separately from feasibility.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
originalPath = path;
testCase.TestData.pathCleanup = onCleanup(@() path(originalPath));
addpath(fileparts(fileparts(mfilename('fullpath'))));
end

function teardownOnce(testCase)
testCase.TestData.pathCleanup = [];
end

function testScaledResidualAndFeasibilityAreIndependent(testCase)
% Two residuals have RMS sqrt(5)*1e-7 and geometric scale 25.
directions = ones(2,1,2);
references = zeros(2,3,2);
references(:,:,1) = [3 4 0;3 4 0];
references(:,:,2) = references(:,:,1);
gij = [1e-7;3e-7];
target = zeros(3,1);
options = struct('ConstraintTolerance',1e-8);
[beta,info] = MiuraPerturbCorrect( ...
    target,directions,references,gij,options);

verifyEqual(testCase,beta,target);
verifyEqual(testCase,info.exitflag,1);
verifyEqual(testCase,info.iterations,1);
verifyEqual(testCase,info.constraintMeasure,sqrt(5)*1e-7/25,'RelTol',1e-12);
verifyEqual(testCase,info.gradientMeasure,4e-6,'RelTol',1e-12);
verifyEqual(testCase,info.objective,5e-14,'RelTol',1e-12);
verifyEqual(testCase,info.relativeStep,Inf);
verifyFalse(testCase,info.constraintSatisfied);

% Tighten only the scaled stopping tolerance while loosening raw feasibility.
options.ResidualTolerance = 1e-10;
options.GradientTolerance = 1e-5;
options.ConstraintTolerance = 1e-6;
[beta,info] = MiuraPerturbCorrect( ...
    target,directions,references,gij,options);
verifyEqual(testCase,beta,target);
verifyEqual(testCase,info.exitflag,2);
verifyGreaterThan(testCase,info.constraintMeasure,options.ResidualTolerance);
verifyTrue(testCase,info.constraintSatisfied);
end

function testScaledGradientUsesTwoNormAndObjective(testCase)
% c=3e-7*x+4e-7*y-1e6: ||gradient||=0.5 and sqrt(2*f)=1e6.
directions = zeros(1,1,2);
directions(:,:,1) = 1e-7;
references = zeros(1,3,2);
references(1,:,2) = [3 4 0];
target = zeros(3,1);
[beta,info] = MiuraPerturbCorrect( ...
    target,directions,references,1e6,[]);

verifyEqual(testCase,beta,target);
verifyEqual(testCase,info.exitflag,2);
verifyEqual(testCase,info.iterations,1);
verifyEqual(testCase,info.gradientMeasure,5e-7,'RelTol',1e-12);
verifyEqual(testCase,info.gradientInfinityNorm,0.4,'RelTol',1e-12);
verifyEqual(testCase,info.constraintMeasure,1,'RelTol',1e-12);
verifyEqual(testCase,info.objective,5e11,'RelTol',1e-12);
verifyEqual(testCase,info.relativeStep,Inf);
verifyFalse(testCase,info.constraintSatisfied);
end

function testRelativeStepKeepsExitFlagThree(testCase)
% c=1e6*x-(2e6-1) at x=2 has residual 1 but relative Newton step about 5e-7.
directions = zeros(1,1,2);
directions(:,:,1) = 1;
references = zeros(1,3,2);
references(1,1,2) = 1e6;
target = [2;0;0];
[beta,info] = MiuraPerturbCorrect( ...
    target,directions,references,2e6-1,[]);

verifyEqual(testCase,info.exitflag,3);
verifyEqual(testCase,info.iterations,1);
verifyEqual(testCase,beta,target);
verifyGreaterThan(testCase,info.constraintMeasure,1e-8);
verifyGreaterThan(testCase,info.gradientMeasure,1e-6);
verifyGreaterThan(testCase,info.relativeStep,0);
verifyLessThanOrEqual(testCase,info.relativeStep,1e-6);
expectedStep = 1e6/(1e12+1e-6+info.finalDamping)/norm(target);
verifyEqual(testCase,info.relativeStep,expectedStep,'RelTol',1e-12);
verifyEqual(testCase,info.objective,0.5,'RelTol',1e-12);
verifyFalse(testCase,info.constraintSatisfied);
end

function testScaleHasUnitFloor(testCase)
% A constant nonzero constraint has zero gradient but is still infeasible.
directions = zeros(1,1,2);
references = zeros(1,3,2);
[~,info] = MiuraPerturbCorrect( ...
    zeros(3,1),directions,references,0.25,[]);
verifyEqual(testCase,info.exitflag,2);
verifyEqual(testCase,info.constraintMeasure,0.25);
verifyEqual(testCase,info.gradientMeasure,0);
verifyFalse(testCase,info.constraintSatisfied);
end

function testIterationLimitReportsFinalAcceptedState(testCase)
% One damped step on c=x-1 lowers f but does not meet these strict thresholds.
directions = zeros(1,1,2);
directions(:,:,1) = 1;
references = zeros(1,3,2);
references(1,1,2) = 1;
options = struct('Regularization',0.25,'MaxIterations',1, ...
    'ResidualTolerance',1e-14,'GradientTolerance',1e-14, ...
    'StepTolerance',1e-14,'ConstraintTolerance',1e-8);
[beta,info] = MiuraPerturbCorrect( ...
    zeros(3,1),directions,references,1,options);
residual = beta(1)-1;
objective = 0.5*residual^2+0.5*options.Regularization*sum(beta.^2);
gradient = [residual;0;0]+options.Regularization*beta;

verifyEqual(testCase,info.exitflag,0);
verifyEqual(testCase,info.iterations,1);
verifyGreaterThan(testCase,beta(1),0);
verifyLessThan(testCase,objective,0.5);
verifyEqual(testCase,info.objective,objective,'RelTol',1e-12);
verifyEqual(testCase,info.constraintMeasure,abs(residual),'RelTol',1e-12);
verifyEqual(testCase,info.gradientMeasure, ...
    norm(gradient)/max(1,sqrt(2*objective)),'AbsTol',1e-14);
verifyGreaterThan(testCase,info.constraintMeasure,options.ResidualTolerance);
verifyGreaterThan(testCase,info.gradientMeasure,options.GradientTolerance);
verifyFalse(testCase,info.constraintSatisfied);
end

function testNonfiniteModelCannotReportConvergence(testCase)
% Finite inputs can still overflow the quadratic model during evaluation.
target = [1e200;0;0];
[beta,info] = MiuraPerturbCorrect( ...
    target,ones(1,1,2),zeros(1,3,2),1,[]);
verifyEqual(testCase,beta,target);
verifyEqual(testCase,info.exitflag,-2);
verifyEqual(testCase,info.iterations,1);
verifyFalse(testCase,isfinite(info.objective));
verifyFalse(testCase,info.constraintSatisfied);
end

function testSparseCurvatureStepAndUpdatedDiagnostics(testCase)
% A coupled two-master problem exercises off-diagonal curvature and a changing J.
directions = zeros(3,2,2);
directions(:,:,1) = [1 2;0.5 -1;1.5 0.4];
directions(:,:,2) = directions(:,:,1).*[1.1 0.9];
references = zeros(3,3,2);
gij = -[1;2;3];
target = [0.2;-0.1;0.05;0.3;-0.15;0.1];
options = struct('Regularization',0.25,'MaxIterations',1, ...
    'ResidualTolerance',1e-14,'GradientTolerance',1e-14,'StepTolerance',1e-14);
[initialResidual,initialJacobian,curvature] = ...
    AffineMetricResidualJacobian(target,directions,references,gij);
% The public evaluator independently assembles curvature by matrix products.
hessian = initialJacobian.'*initialJacobian+curvature+options.Regularization*eye(6);
hessian = 0.5*(hessian+hessian.');
hessianScale = max(1,norm(hessian,inf));
mu = max([1e-8*hessianScale,-min(eig(hessian))+sqrt(eps)*hessianScale,eps*hessianScale]);
expected = target-(hessian+mu*eye(6))\(initialJacobian.'*initialResidual);
[beta,info] = MiuraPerturbCorrect(target,directions,references,gij,options);
verifyEqual(testCase,info.exitflag,0);
verifyEqual(testCase,beta,expected,'AbsTol',2e-13);
[residual,jacobian] = AffineMetricResidualJacobian(beta,directions,references,gij);
gradient = jacobian.'*residual+options.Regularization*(beta-target);
objective = 0.5*sum(residual.^2)+0.5*options.Regularization*sum((beta-target).^2);
verifyLessThan(testCase,objective,0.5*sum(initialResidual.^2));
verifyGreaterThan(testCase,norm(jacobian-initialJacobian,'fro'),0.1);
verifyEqual(testCase,info.objective,objective,'AbsTol',1e-13);
verifyEqual(testCase,info.maxResidual,norm(residual,inf),'AbsTol',1e-13);
verifyEqual(testCase,info.gradientInfinityNorm,norm(gradient,inf),'AbsTol',1e-13);
end

function testCurvatureCacheTracksModelChanges(testCase)
% Reusing a cache must match a fresh solve after values or dimensions change.
directions = zeros(2,2,2);
directions(:,:,1) = [1 0.3;0.2 1];
directions(:,:,2) = [0.7 0.5;0.4 1.2];
references = zeros(2,3,2);
target = [0.2;-0.1;0.05;0.3;-0.15;0.1];
options = struct('Regularization',0.25,'MaxIterations',1, ...
    'ResidualTolerance',1e-14,'GradientTolerance',1e-14,'StepTolerance',1e-14);
for variant = 1:4
    if variant == 2
        directions(1,1,2) = 1.4;
    elseif variant == 3
        directions = directions(:,1,:);
        target = target(1:3);
    elseif variant == 4
        directions(:,:,2) = 0;
        references(:,:,2) = [1 0.5 0;0.2 1 0.4];
    end
    [cachedBeta,cachedInfo] = MiuraPerturbCorrect(target,directions,references,-[1;2],options);
    clear MiuraPerturbCorrect
    [freshBeta,freshInfo] = MiuraPerturbCorrect(target,directions,references,-[1;2],options);
    verifyEqual(testCase,cachedBeta,freshBeta);
    verifyEqual(testCase,cachedInfo,freshInfo);
end
end

function testRejectedTrialDoesNotOverwriteAcceptedState(testCase)
% Negative curvature produces an excessive first trial; allow only that trial.
target = [0.1;0;0];
directions = ones(1,1,2);
references = zeros(1,3,2);
options = struct('MaxInnerIterations',1);
[beta,info] = MiuraPerturbCorrect(target,directions,references,1,options);
[residual,jacobian] = AffineMetricResidualJacobian(target,directions,references,1);
verifyEqual(testCase,info.exitflag,-1);
verifyEqual(testCase,info.iterations,1);
verifyEqual(testCase,beta,target);
verifyEqual(testCase,info.maxResidual,norm(residual,inf),'AbsTol',1e-14);
verifyEqual(testCase,info.objective,0.5*sum(residual.^2),'AbsTol',1e-14);
verifyEqual(testCase,info.gradientInfinityNorm,norm(jacobian.'*residual,inf),'AbsTol',1e-14);
end
