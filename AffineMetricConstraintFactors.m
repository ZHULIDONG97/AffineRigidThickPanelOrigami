function [constraintDirections,referenceComponents,gij, ...
    selectedConstraints,candidateConstraints,recoveryMatrices] = ...
    AffineMetricConstraintFactors(B,x0,y0,z0,panelNodes, ...
    intrinsicDimensions,tolerance)
%AFFINEMETRICCONSTRAINTFACTORS Build and preselect A'*A=I constraints.

if nargin < 7 || isempty(tolerance)
    tolerance = 1e-10;
end

% Validate the shared affine coordinates before constructing panel maps.
x0 = x0(:);
y0 = y0(:);
z0 = z0(:);
nv = numel(x0);
if numel(y0) ~= nv || numel(z0) ~= nv || size(B,1) ~= nv
    error('AffineMetricConstraintFactors:CoordinateSizeMismatch', ...
        'B and the three reference-coordinate vectors must use the same nodes.');
end
if isempty(B) || any(~isfinite(B(:))) || ...
        any(~isfinite([x0;y0;z0])) || ~isreal(B)
    error('AffineMetricConstraintFactors:InvalidAffineBasis', ...
        'B must be a finite, real affine basis with at least one column.');
end
if ~iscell(panelNodes) || isempty(panelNodes)
    error('AffineMetricConstraintFactors:InvalidPanels', ...
        'panelNodes must be a nonempty cell array of node-index vectors.');
end
if ~isscalar(tolerance) || ~isfinite(tolerance) || tolerance <= 0
    error('AffineMetricConstraintFactors:InvalidTolerance', ...
        'tolerance must be a finite positive scalar.');
end

nPanels = numel(panelNodes);
if isscalar(intrinsicDimensions)
    intrinsicDimensions = repmat(intrinsicDimensions,nPanels,1);
else
    intrinsicDimensions = intrinsicDimensions(:);
end
if numel(intrinsicDimensions) ~= nPanels || ...
        any(~ismember(intrinsicDimensions,[2 3]))
    error('AffineMetricConstraintFactors:InvalidIntrinsicDimension', ...
        'Use intrinsic dimension 2 or 3 for every panel.');
end

coordinates0 = [x0,y0,z0];
nCandidates = sum(intrinsicDimensions.*(intrinsicDimensions+1)/2);
candidateU = zeros(nCandidates,nv);
candidateV = zeros(nCandidates,nv);
candidateTarget = zeros(nCandidates,1);
candidateConstraints = zeros(nCandidates,3);
recoveryMatrices = cell(nPanels,1);
candidateIndex = 0;

for panelIndex = 1:nPanels
    % Center the reference panel so its recovery map excludes translation.
    nodes = unique(panelNodes{panelIndex}(:),'stable');
    d = intrinsicDimensions(panelIndex);
    if numel(nodes) < d+1 || any(nodes < 1 | nodes > nv) || ...
            any(nodes ~= round(nodes))
        error('AffineMetricConstraintFactors:InvalidPanelNodes', ...
            'Panel %d needs at least %d valid, distinct nodes.', ...
            panelIndex,d+1);
    end
    X0 = coordinates0(nodes,:).';
    H = eye(numel(nodes)) ...
        -ones(numel(nodes),numel(nodes))/numel(nodes);
    Y0 = X0*H;

    % Use the panel plane for d=2 and the centered 3-D frame for d=3.
    [referenceBasis,singularValueMatrix] = svd(Y0,'econ');
    singularValues = diag(singularValueMatrix);
    if isempty(singularValues) || singularValues(1) == 0
        referenceRank = 0;
    else
        rankTolerance = tolerance*max(size(Y0))*singularValues(1);
        referenceRank = sum(singularValues > rankTolerance);
    end
    if (d == 2 && referenceRank ~= 2) || (d == 3 && referenceRank < 3)
        error('AffineMetricConstraintFactors:DegeneratePanel', ...
            ['Panel %d has centered rank %d, but intrinsic dimension ' ...
             '%d was requested.'],panelIndex,referenceRank,d);
    end
    if d == 2
        Xi = referenceBasis(:,1:2).'*Y0;
    else
        Xi = Y0;
    end

    % Recover A=X_e*R_e using only the fixed reference geometry.
    recovery = H*Xi.'/(Xi*Xi.');
    recoveryMatrices{panelIndex} = recovery;
    referenceMap = X0*recovery;
    if norm(referenceMap.'*referenceMap-eye(d),'fro') > 100*tolerance
        error('AffineMetricConstraintFactors:ReferenceMetricMismatch', ...
            'Panel %d does not recover the identity reference metric.', ...
            panelIndex);
    end

    % Ensure A captures every allowed panel motion, not only its projection.
    affineCoordinates = [ones(numel(nodes),1),Xi.'];
    panelBasis = B(nodes,:);
    affineBasisResidual = panelBasis ...
        -affineCoordinates*(affineCoordinates\panelBasis);
    if norm(affineBasisResidual,'fro') ...
            > 100*tolerance*max(1,norm(panelBasis,'fro')) || ...
            norm(Y0-referenceMap*Xi,'fro') ...
            > 100*tolerance*max(1,norm(Y0,'fro'))
        error('AffineMetricConstraintFactors:NonaffinePanelMotion', ...
            ['Panel %d or its supplied basis contains motion outside ' ...
             'the recovered affine panel space.'],panelIndex);
    end

    % Preset every upper-triangular component of A_e'*A_e-I.
    globalRecovery = zeros(nv,d);
    globalRecovery(nodes,:) = recovery;
    for alpha = 1:d
        for beta = alpha:d
            candidateIndex = candidateIndex+1;
            candidateU(candidateIndex,:) = globalRecovery(:,alpha).';
            candidateV(candidateIndex,:) = globalRecovery(:,beta).';
            candidateTarget(candidateIndex) = double(alpha == beta);
            candidateConstraints(candidateIndex,:) = ...
                [panelIndex,alpha,beta];
        end
    end
end

% Project the fixed recovery vectors into the simulator's affine space.
directionU = candidateU*B;
directionV = candidateV*B;
referenceU = candidateU*coordinates0;
referenceV = candidateV*coordinates0;
referenceResidual = sum(referenceU.*referenceV,2)-candidateTarget;

% Include both quadratic coefficients and the reference Jacobian in rank tests.
uu = directionU*directionU.';
vv = directionV*directionV.';
uv = directionU*directionV.';
quadraticGram = 0.5*(uu.*vv+uv.*uv.');
referenceJacobian = [ ...
    referenceV(:,1).*directionU+referenceU(:,1).*directionV, ...
    referenceV(:,2).*directionU+referenceU(:,2).*directionV, ...
    referenceV(:,3).*directionU+referenceU(:,3).*directionV];
polynomialGram = 3*quadraticGram ...
    +referenceJacobian*referenceJacobian.' ...
    +referenceResidual*referenceResidual.';
polynomialGram = 0.5*(polynomialGram+polynomialGram.');

% Normalize the candidate polynomials before rank-revealing selection.
candidateScale = sqrt(max(diag(polynomialGram),0));
largestScale = max(candidateScale);
if largestScale == 0
    error('AffineMetricConstraintFactors:NoActiveConstraint', ...
        'Every preset metric constraint is constant in the supplied affine space.');
end
activeCandidates = find(candidateScale > tolerance*largestScale);
activeScale = candidateScale(activeCandidates);
normalizedGram = polynomialGram(activeCandidates,activeCandidates) ...
    ./(activeScale*activeScale.');
normalizedGram = 0.5*(normalizedGram+normalizedGram.');

% Factor the small Gram matrix, then use CPQR to retain independent rows.
[gramVectors,gramValues] = eig(normalizedGram,'vector');
[gramValues,gramOrder] = sort(real(gramValues),'descend');
gramVectors = real(gramVectors(:,gramOrder));
largestEigenvalue = max(gramValues(1),0);
rankTolerance = max(tolerance^2*largestEigenvalue, ...
    numel(activeCandidates)*eps(largestEigenvalue));
constraintRank = sum(gramValues > rankTolerance);
if constraintRank == 0
    error('AffineMetricConstraintFactors:NoIndependentConstraint', ...
        'No independent metric constraint remains after rank screening.');
end
coefficientFactor = diag(sqrt(gramValues(1:constraintRank))) ...
    *gramVectors(:,1:constraintRank).';
[~,~,constraintOrder] = qr(coefficientFactor,'vector');
selectedCandidateIndices = sort(activeCandidates( ...
    constraintOrder(1:constraintRank)));

% Cache only the selected low-rank factors for every later solver call.
selectedConstraints = candidateConstraints(selectedCandidateIndices,:);
constraintDirections = zeros(constraintRank,size(B,2),2);
constraintDirections(:,:,1) = directionU(selectedCandidateIndices,:);
constraintDirections(:,:,2) = directionV(selectedCandidateIndices,:);
referenceComponents = zeros(constraintRank,3,2);
referenceComponents(:,:,1) = referenceU(selectedCandidateIndices,:);
referenceComponents(:,:,2) = referenceV(selectedCandidateIndices,:);
gij = -referenceResidual(selectedCandidateIndices);
end
