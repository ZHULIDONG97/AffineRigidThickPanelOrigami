function tests = test_low_rank_derivatives
% Check metric derivatives against independent scalar finite differences.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
originalPath = path;
originalRng = rng;
testCase.TestData.pathCleanup = onCleanup(@() path(originalPath));
testCase.TestData.rngCleanup = onCleanup(@() rng(originalRng));
addpath(fileparts(fileparts(mfilename('fullpath'))));
rng(20260919,'twister');
testCase.TestData.directions = 0.5*randn(7,3,2);
testCase.TestData.references = randn(7,3,2);
testCase.TestData.gij = 0.3*randn(7,1);
testCase.TestData.beta = 0.2*randn(9,1);
end

function teardownOnce(testCase)
% Release the cleanup objects even when an individual verification fails.
testCase.TestData.pathCleanup = [];
testCase.TestData.rngCleanup = [];
end

function testResidualMatchesScalarGeometry(testCase)
data = testCase.TestData;
actual = AffineMetricResidualJacobian( ...
    data.beta,data.directions,data.references,data.gij);
expected = scalarResidual( ...
    data.beta,data.directions,data.references,data.gij);
verifyEqual(testCase,actual,expected,'AbsTol',5e-14);

% The public function accepts either orientation of the coordinate vector.
rowInput = AffineMetricResidualJacobian( ...
    data.beta.',data.directions,data.references,data.gij.');
verifyEqual(testCase,rowInput,actual,'AbsTol',5e-14);
end

function testPenaltyUsesOnlyMetricResidual(testCase)
data = testCase.TestData;
actual = equalityPenaltyObjective(data.beta,1e6,data.directions,data.references,data.gij);
expected = 0.5e6*sum(scalarResidual(data.beta,data.directions,data.references,data.gij).^2);
verifyEqual(testCase,actual,expected,'RelTol',1e-13);
end

function testJacobianMatchesCentralDifference(testCase)
data = testCase.TestData;
[~,jacobian] = AffineMetricResidualJacobian( ...
    data.beta,data.directions,data.references,data.gij);
step = 1e-5;
finiteDifference = zeros(size(jacobian));
for dof = 1:numel(data.beta)
    offset = zeros(size(data.beta));
    offset(dof) = step;
    plusResidual = scalarResidual( ...
        data.beta+offset,data.directions,data.references,data.gij);
    minusResidual = scalarResidual( ...
        data.beta-offset,data.directions,data.references,data.gij);
    finiteDifference(:,dof) = (plusResidual-minusResidual)/(2*step);
end
relativeError = norm(jacobian-finiteDifference,'fro') ...
    /max(1,norm(finiteDifference,'fro'));
verifyLessThan(testCase,relativeError,2e-9);
end

function testExactObjectiveHessianMatchesCentralDifference(testCase)
data = testCase.TestData;
[~,jacobian,residualCurvature] = AffineMetricResidualJacobian( ...
    data.beta,data.directions,data.references,data.gij);
hessian = jacobian.'*jacobian+residualCurvature;
objective = @(beta) 0.5*sum(scalarResidual( ...
    beta,data.directions,data.references,data.gij).^2);
step = 1e-4;
nDof = numel(data.beta);
finiteDifference = zeros(nDof,nDof);
initialObjective = objective(data.beta);
for firstDof = 1:nDof
    firstOffset = zeros(nDof,1);
    firstOffset(firstDof) = step;
    finiteDifference(firstDof,firstDof) = ...
        (objective(data.beta+firstOffset)-2*initialObjective ...
        +objective(data.beta-firstOffset))/step^2;
    for secondDof = firstDof+1:nDof
        secondOffset = zeros(nDof,1);
        secondOffset(secondDof) = step;
        entry = (objective(data.beta+firstOffset+secondOffset) ...
            -objective(data.beta+firstOffset-secondOffset) ...
            -objective(data.beta-firstOffset+secondOffset) ...
            +objective(data.beta-firstOffset-secondOffset))/(4*step^2);
        finiteDifference(firstDof,secondDof) = entry;
        finiteDifference(secondDof,firstDof) = entry;
    end
end
relativeError = norm(hessian-finiteDifference,'fro') ...
    /max(1,norm(finiteDifference,'fro'));
verifyLessThan(testCase,relativeError,2e-6);
verifyEqual(testCase,hessian,hessian.','AbsTol',1e-13);
% Ensure the fixture exercises the exact curvature beyond Gauss-Newton.
verifyGreaterThan(testCase,norm(residualCurvature,'fro'),0.1);
end

function testRejectsInconsistentDimensions(testCase)
data = testCase.TestData;
errorId = 'AffineMetricResidualJacobian:InvalidConstraintData';
verifyError(testCase,@() AffineMetricResidualJacobian( ...
    data.beta(1:end-1),data.directions,data.references,data.gij),errorId);
verifyError(testCase,@() AffineMetricResidualJacobian( ...
    data.beta,data.directions(:,:,1),data.references,data.gij),errorId);
verifyError(testCase,@() AffineMetricResidualJacobian( ...
    data.beta,data.directions,data.references(:,1:2,:),data.gij),errorId);
verifyError(testCase,@() AffineMetricResidualJacobian( ...
    data.beta,data.directions,data.references,data.gij(1:end-1)),errorId);
end

function testRejectsNonfiniteInputs(testCase)
data = testCase.TestData;
for invalidValue = [NaN Inf]
    for inputIndex = 1:4
        inputs = {data.beta,data.directions,data.references,data.gij};
        inputs{inputIndex}(1) = invalidValue;
        verifyError(testCase,@() AffineMetricResidualJacobian(inputs{:}), ...
            'AffineMetricResidualJacobian:NonfiniteInput');
    end
end
end

function residual = scalarResidual(beta,directions,references,gij)
% Form current and reference dot products without the production evaluator.
nConstraints = size(directions,1);
nMasters = size(directions,2);
residual = -gij(:);
for constraint = 1:nConstraints
    for coordinate = 1:3
        indices = (coordinate-1)*nMasters+(1:nMasters);
        referenceU = references(constraint,coordinate,1);
        referenceV = references(constraint,coordinate,2);
        currentU = referenceU+directions(constraint,:,1)*beta(indices);
        currentV = referenceV+directions(constraint,:,2)*beta(indices);
        residual(constraint) = residual(constraint) ...
            +currentU*currentV-referenceU*referenceV;
    end
end
end
