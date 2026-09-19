% Generate a relative 1-degree target for node 18, then correct its perturbation once by LM.
clear;
clc;
close all;

codeDir = fileparts(mfilename('fullpath'));
modelFile = fullfile(codeDir,'four_panel_reference.mat');

%% Settings
dragNode = 18;
dragRotationEdge = [17 24];
dragAngleStep = deg2rad(1);
regularization = 1e-7;
maximumRotationAngle = deg2rad(180);
angleCompletionTolerance = deg2rad(0.001);
constraintRankTolerance = 1e-10;

if ~isscalar(dragAngleStep) || ~isfinite(dragAngleStep) || dragAngleStep <= 0 || ...
        regularization <= 0 || ...
        maximumRotationAngle <= 0 || ~isfinite(maximumRotationAngle) || ...
        ~isfinite(angleCompletionTolerance) || angleCompletionTolerance <= 0 || ...
        angleCompletionTolerance >= maximumRotationAngle
    error('demo_miura_drag:InvalidDrag', ...
        ['Angular increment and regularization must be positive; ' ...
         'the finite angle tolerance must lie between zero and the angle limit.']);
end

%% Zero-displacement model
model = load(modelFile,'x0','y0','z0','B','T', ...
    'mainPointIndices','fix','panels');
x0 = model.x0(:);
y0 = model.y0(:);
z0 = model.z0(:);
mainPointIndices = model.mainPointIndices(:).';
fixedNodes = model.fix(:).';
panels = model.panels;

mappedBFull = model.T\model.B;
[~,fixedMasterLocalIndices] = ismember(fixedNodes,mainPointIndices);
freeMasterLocalIndices = setdiff(1:numel(mainPointIndices), ...
    fixedMasterLocalIndices,'stable');
dragMasterFreeIndex = find( ...
    mainPointIndices(freeMasterLocalIndices) == dragNode,1);
if isempty(dragMasterFreeIndex)
    error('demo_miura_drag:InvalidDragNode', ...
        'dragNode must be an unfixed master node.');
end

mappedB = mappedBFull(:,freeMasterLocalIndices);
nFreeMasters = size(mappedB,2);
dragDofIndices = dragMasterFreeIndex+[0,nFreeMasters,2*nFreeMasters];

% Select panel metric constraints once in the full compatible affine space.
panelNodes = cellfun(@(panelFaces) unique([panelFaces{:}],'stable'), ...
    panels,'UniformOutput',false);
[rigidityDirectionsFull,rigidityReferences,rigidityOffsets, ...
    metricConstraints,candidateMetricConstraints] = ...
    AffineMetricConstraintFactors(mappedBFull,x0,y0,z0,panelNodes,3, ...
    constraintRankTolerance);
rigidityDirections = ...
    rigidityDirectionsFull(:,freeMasterLocalIndices,:);
fprintf('Thick-panel metric constraints: %d preset, %d independent.\n', ...
    size(candidateMetricConstraints,1),size(metricConstraints,1));

% Nodes 17-18-24 form the broad surface containing the dragged node.
normalFaceNodes = [17 18 24];

% Each row defines the two panels adjacent to one physical hinge.
rotationNodeIndices = [16 17 24 18; ...
    11 12 21 13; ...
    14 15 20 10; ...
    2 19 1 9];
rotationEdges = rotationNodeIndices(:,2:3);
rotationAngleTolerance = deg2rad(1e-8);
angleBoundaryTolerance = deg2rad(0.1);
maximumPerturbationBisections = 30;

%% Perturbation correction
nNodes = numel(x0);
if numel(dragRotationEdge) ~= 2 || numel(unique(dragRotationEdge)) ~= 2 || ...
        any(dragRotationEdge < 1 | dragRotationEdge > nNodes ...
        | dragRotationEdge ~= round(dragRotationEdge)) || ismember(dragNode,dragRotationEdge)
    error('demo_miura_drag:InvalidDragAxis', ...
        'The drag rotation axis needs two distinct valid nodes excluding the dragged node.');
end
if any(rotationNodeIndices(:) < 1 | rotationNodeIndices(:) > nNodes)
    error('demo_miura_drag:InvalidRotationNodes', ...
        'A rotation-angle constraint contains an invalid node index.');
end
rotationAngles = RotationAngle( ...
    size(rotationNodeIndices,1),rotationNodeIndices,x0,y0,z0);
if any(~isfinite(rotationAngles)) || any(abs(rotationAngles) > ...
        maximumRotationAngle+rotationAngleTolerance)
    error('demo_miura_drag:InitialRotationViolation', ...
        'The initial model already exceeds the rotation-angle limit.');
end

% Grow storage as needed; its initial capacity does not limit the step count.
historyCapacity = 64;
xHistory = zeros(nNodes,historyCapacity+1);
yHistory = zeros(nNodes,historyCapacity+1);
zHistory = zeros(nNodes,historyCapacity+1);
xHistory(:,1) = x0;
yHistory(:,1) = y0;
zHistory(:,1) = z0;

coordinates = [x0,y0,z0];
masterDisplacement = zeros(3*nFreeMasters,1);
previousPanelNormal = [];
maximumResidual = 0;
maximumAcceptedRotationAngle = max(abs(rotationAngles));
angleLimitReached = maximumRotationAngle-max(abs(rotationAngles)) <= angleCompletionTolerance;
attemptedSteps = 0;
completedSteps = 0;
rejectedRotationAngle = NaN;
rejectedRotationEdge = [NaN NaN];
solverOptions = struct('Regularization',regularization, ...
    'ConstraintTolerance',1e-8,'MaxIterations',150);

% Time the complete drag loop; rejected LM trials count in the step average.
stepWallTime = zeros(historyCapacity,1);
totalLmCalls = 0;
totalLmIterations = 0;
failedCorrections = 0;
angleLimitTrials = 0;
perturbationBisections = 0;
perturbationFraction = 1;
perturbationScale = 1;
constraintToleranceExceedances = 0;
solveCpuStart = cputime;
solveTimer = tic;
while ~angleLimitReached
    stepIndex = attemptedSteps+1;
    if stepIndex > historyCapacity
        historyCapacity = 2*historyCapacity;
        xHistory(:,historyCapacity+1) = 0;
        yHistory(:,historyCapacity+1) = 0;
        zHistory(:,historyCapacity+1) = 0;
        stepWallTime(historyCapacity,1) = 0;
    end
    stepTimer = tic;
    faceCoordinates = coordinates(normalFaceNodes,:);
    panelNormal = cross( ...
        faceCoordinates(2,:)-faceCoordinates(1,:), ...
        faceCoordinates(3,:)-faceCoordinates(1,:));
    panelNormal = panelNormal/norm(panelNormal);
    if ~isempty(previousPanelNormal) && ...
            dot(panelNormal,previousPanelNormal) < 0
        panelNormal = -panelNormal;
    end

    % Rotate around the current hinge, choosing the old outward drag direction.
    hingeOrigin = coordinates(dragRotationEdge(1),:);
    hingeAxis = coordinates(dragRotationEdge(2),:)-hingeOrigin;
    hingeAxisLength = norm(hingeAxis);
    if ~isfinite(hingeAxisLength) || hingeAxisLength == 0
        error('demo_miura_drag:DegenerateDragAxis','The drag hinge has zero or invalid length.');
    end
    hingeAxis = hingeAxis/hingeAxisLength;
    radialPosition = coordinates(dragNode,:)-hingeOrigin;
    if dot(cross(hingeAxis,radialPosition),panelNormal) < 0
        hingeAxis = -hingeAxis;
    end

    % Generate a target 1 degree beyond the current corrected state, not stepIndex degrees.
    % Convert that relative rotation into a Cartesian perturbation of the drag node.
    targetDragCoordinates = hingeOrigin+radialPosition*cos(dragAngleStep) ...
        +cross(hingeAxis,radialPosition)*sin(dragAngleStep) ...
        +hingeAxis*dot(hingeAxis,radialPosition)*(1-cos(dragAngleStep));
    dragCoordinatePerturbation = targetDragCoordinates-coordinates(dragNode,:);
    % Bisect only the predicted coordinates. These angle checks do not call LM.
    lowerFraction = 0;
    upperFraction = perturbationScale;
    perturbationFraction = perturbationScale;
    targetMasterDisplacement = masterDisplacement;
    for predictionIndex = 0:maximumPerturbationBisections
        predictedMasterDisplacement = masterDisplacement;
        predictedMasterDisplacement(dragDofIndices) = ...
            predictedMasterDisplacement(dragDofIndices) ...
            +perturbationFraction*dragCoordinatePerturbation.';
        predictedCoordinates = [x0,y0,z0] ...
            +mappedB*reshape(predictedMasterDisplacement,nFreeMasters,3);
        predictedPrincipalAngles = RotationAngle( ...
            size(rotationNodeIndices,1),rotationNodeIndices, ...
            predictedCoordinates(:,1),predictedCoordinates(:,2),predictedCoordinates(:,3));
        predictedRotationAngles = rotationAngles+atan2( ...
            sin(predictedPrincipalAngles-rotationAngles), ...
            cos(predictedPrincipalAngles-rotationAngles));
        predictedMaximumAngle = max(abs(predictedRotationAngles));
        if all(isfinite(predictedRotationAngles)) && ...
                predictedMaximumAngle <= maximumRotationAngle+rotationAngleTolerance
            lowerFraction = perturbationFraction;
            targetMasterDisplacement = predictedMasterDisplacement;
            if perturbationFraction == perturbationScale || ...
                    maximumRotationAngle-predictedMaximumAngle <= angleBoundaryTolerance
                break
            end
        else
            upperFraction = perturbationFraction;
        end
        if predictionIndex < maximumPerturbationBisections
            perturbationFraction = 0.5*(lowerFraction+upperFraction);
            perturbationBisections = perturbationBisections+1;
        end
    end
    perturbationFraction = lowerFraction;
    if perturbationFraction == 0
        fprintf(['Step %d stopped: no admissible positive perturbation ' ...
            'was found within the angle-search limit.\n'],stepIndex);
        break
    end
    [trialMasterDisplacement,solverInfo] = MiuraPerturbCorrect( ...
        targetMasterDisplacement,rigidityDirections, ...
        rigidityReferences,rigidityOffsets,solverOptions);
    attemptedSteps = stepIndex;
    totalLmCalls = totalLmCalls+1;
    totalLmIterations = totalLmIterations+solverInfo.iterations;

    % Match the reference workflow: stop on LM failure, without another solve.
    if solverInfo.exitflag <= 0 || ...
            any(~isfinite(trialMasterDisplacement)) || ~isfinite(solverInfo.maxResidual)
        failedCorrections = failedCorrections+1;
        fprintf(['Step %d stopped after one LM call: exitflag=%d, ' ...
            'max|c|=%.3e. The previous state was retained.\n'], ...
            stepIndex,solverInfo.exitflag,solverInfo.maxResidual);
        stepWallTime(stepIndex) = toc(stepTimer);
        break
    end

    trialCoordinates = [x0,y0,z0] ...
        +mappedB*reshape(trialMasterDisplacement,nFreeMasters,3);
    principalRotationAngles = RotationAngle( ...
        size(rotationNodeIndices,1),rotationNodeIndices, ...
        trialCoordinates(:,1),trialCoordinates(:,2),trialCoordinates(:,3));
    rotationAngleIncrement = atan2( ...
        sin(principalRotationAngles-rotationAngles), ...
        cos(principalRotationAngles-rotationAngles));
    trialRotationAngles = rotationAngles+rotationAngleIncrement;
    if any(~isfinite(trialCoordinates(:))) || any(~isfinite(trialRotationAngles))
        failedCorrections = failedCorrections+1;
        fprintf('Step %d stopped: corrected geometry is nonfinite; the previous state was retained.\n',stepIndex);
        stepWallTime(stepIndex) = toc(stepTimer);
        break
    end
    [trialMaximumRotationAngle,rotationIndex] = max(abs(trialRotationAngles));

    % Keep an overshoot out of the history and halve the next step's coordinates.
    if trialMaximumRotationAngle > maximumRotationAngle+rotationAngleTolerance
        angleLimitTrials = angleLimitTrials+1;
        rejectedRotationAngle = trialMaximumRotationAngle;
        rejectedRotationEdge = rotationEdges(rotationIndex,:);
        fprintf(['Step %d rejected at hinge [%d %d]: trial %.6f deg exceeds ' ...
            'the %.6f deg limit. Halving the next coordinate perturbation.\n'], ...
            stepIndex,rejectedRotationEdge(1),rejectedRotationEdge(2), ...
            rad2deg(rejectedRotationAngle),rad2deg(maximumRotationAngle));
        stepWallTime(stepIndex) = toc(stepTimer);
        perturbationScale = 0.5*perturbationFraction;
        continue
    end

    % An unchanged state would repeat the same target and deterministic LM solve.
    angleLimitReached = maximumRotationAngle-trialMaximumRotationAngle <= angleCompletionTolerance;
    if ~angleLimitReached && isequal(trialMasterDisplacement,masterDisplacement)
        failedCorrections = failedCorrections+1;
        fprintf('Step %d stopped: LM returned an unchanged state before the angle target.\n',stepIndex);
        stepWallTime(stepIndex) = toc(stepTimer);
        break
    end

    % Accept positive LM termination; the raw residual remains a diagnostic.
    masterDisplacement = trialMasterDisplacement;
    coordinates = trialCoordinates;
    previousPanelNormal = panelNormal;
    rotationAngles = trialRotationAngles;
    completedSteps = completedSteps+1;
    xHistory(:,completedSteps+1) = coordinates(:,1);
    yHistory(:,completedSteps+1) = coordinates(:,2);
    zHistory(:,completedSteps+1) = coordinates(:,3);
    maximumResidual = max(maximumResidual,solverInfo.maxResidual);
    maximumAcceptedRotationAngle = max( ...
        maximumAcceptedRotationAngle,trialMaximumRotationAngle);
    constraintToleranceExceedances = constraintToleranceExceedances ...
        +double(~solverInfo.constraintSatisfied);
    fprintf(['Step %d: coordinate perturbation fraction %.6g, one LM call, %d iterations, ' ...
        'max|c| %.3e, max angle %.6f deg\n'], ...
        stepIndex,perturbationFraction, ...
        solverInfo.iterations,solverInfo.maxResidual, ...
        rad2deg(trialMaximumRotationAngle));
    stepWallTime(stepIndex) = toc(stepTimer);
end
solveTime = toc(solveTimer);
solveCpuTime = cputime-solveCpuStart;
stepWallTime = stepWallTime(1:attemptedSteps);
averageStepWallTime = NaN;
averageStepCpuTime = NaN;
if attemptedSteps > 0
    averageStepWallTime = solveTime/attemptedSteps;
    averageStepCpuTime = solveCpuTime/attemptedSteps;
end

xHistory = xHistory(:,1:completedSteps+1);
yHistory = yHistory(:,1:completedSteps+1);
zHistory = zHistory(:,1:completedSteps+1);
x_fs = xHistory(:,end);
y_fs = yHistory(:,end);
z_fs = zHistory(:,end);
fprintf(['Finished %d accepted / %d attempted steps in %.6f s simulation time; ' ...
    'max|c| %.3e, max accepted angle %.6f deg.\n'], ...
    completedSteps,attemptedSteps,solveTime, ...
    maximumResidual,rad2deg(maximumAcceptedRotationAngle));
if angleLimitReached
    fprintf('Reached %.6f deg within %.6f deg of the angle target.\n', ...
        rad2deg(max(abs(rotationAngles))),rad2deg(angleCompletionTolerance));
else
    fprintf('The angle target was not reached; the last accepted state is retained.\n');
end

% Includes prediction, full LM calls, geometry checks, history, and step output.
fprintf('Simulation loop including outer work: %.6f s wall / %.6f s CPU.\n',solveTime,solveCpuTime);
if attemptedSteps > 0
    fprintf(['Mean per attempted step (rejections included): %.3f ms wall / %.3f ms CPU; ' ...
        'wall median %.3f ms, range %.3f-%.3f ms.\n'], ...
        1000*averageStepWallTime,1000*averageStepCpuTime, ...
        1000*median(stepWallTime),1000*min(stepWallTime),1000*max(stepWallTime));
    fprintf('Mean LM outer iterations per call: %.2f.\n',totalLmIterations/totalLmCalls);
end
fprintf(['Solver work: %d LM calls; %d LM outer iterations; %d failed corrections; ' ...
    '%d rejected angle-limit trials; %d perturbation-coordinate bisections.\n'], ...
    totalLmCalls,totalLmIterations,failedCorrections,angleLimitTrials,perturbationBisections);
fprintf(['Accepted steps with raw max|c| above %.3e: %d. ' ...
    'Residuals are diagnostic and do not trigger retries.\n'], ...
    solverOptions.ConstraintTolerance,constraintToleranceExceedances);

%% Animation
xLimits = [min(xHistory(:)),max(xHistory(:))];
yLimits = [min(yHistory(:)),max(yHistory(:))];
zLimits = [min(zHistory(:)),max(zHistory(:))];
padding = 0.1*max([diff(xLimits),diff(yLimits),diff(zLimits)]);
xLimits = xLimits+[-padding padding];
yLimits = yLimits+[-padding padding];
zLimits = zLimits+[-padding padding];

animationFigure = figure('Color','w','Name','Miura node 18 drag', ...
    'NumberTitle','off');
animationAxes = axes(animationFigure);
axes(animationAxes);
panelHandles = DrawSolidPanels(xHistory(:,1),yHistory(:,1), ...
    zHistory(:,1),panels,'FaceAlpha',0.56,'LineWidth',1.0);
faceNodes = [panels{:}];
hold(animationAxes,'on');
fixedHandle = scatter3(animationAxes,xHistory(fixedNodes,1), ...
    yHistory(fixedNodes,1),zHistory(fixedNodes,1),80, ...
    [0.08 0.08 0.08],'filled','MarkerEdgeColor','w');
dragHandle = scatter3(animationAxes,xHistory(dragNode,1), ...
    yHistory(dragNode,1),zHistory(dragNode,1),150, ...
    [0.02 0.48 0.95],'filled','MarkerEdgeColor','w');
faceCoordinates = [xHistory(normalFaceNodes,1), ...
    yHistory(normalFaceNodes,1),zHistory(normalFaceNodes,1)];
arrowDirection = cross( ...
    faceCoordinates(2,:)-faceCoordinates(1,:), ...
    faceCoordinates(3,:)-faceCoordinates(1,:));
arrowDirection = padding*arrowDirection/norm(arrowDirection);
tractionHandle = quiver3(animationAxes,xHistory(dragNode,1), ...
    yHistory(dragNode,1),zHistory(dragNode,1), ...
    arrowDirection(1),arrowDirection(2),arrowDirection(3),0, ...
    'Color','r','LineWidth',2.5,'MaxHeadSize',0.8);
hold(animationAxes,'off');
axis(animationAxes,'equal');
axis(animationAxes,'vis3d');
axis(animationAxes,'off');
view(animationAxes,210,24);
xlim(animationAxes,xLimits);
ylim(animationAxes,yLimits);
zlim(animationAxes,zLimits);

for historyIndex = 1:completedSteps+1
    for faceIndex = 1:numel(faceNodes)
        nodes = faceNodes{faceIndex};
        set(panelHandles(faceIndex), ...
            'XData',xHistory(nodes,historyIndex), ...
            'YData',yHistory(nodes,historyIndex), ...
            'ZData',zHistory(nodes,historyIndex));
    end
    set(fixedHandle, ...
        'XData',xHistory(fixedNodes,historyIndex), ...
        'YData',yHistory(fixedNodes,historyIndex), ...
        'ZData',zHistory(fixedNodes,historyIndex));
    set(dragHandle, ...
        'XData',xHistory(dragNode,historyIndex), ...
        'YData',yHistory(dragNode,historyIndex), ...
        'ZData',zHistory(dragNode,historyIndex));
    faceCoordinates = [xHistory(normalFaceNodes,historyIndex), ...
        yHistory(normalFaceNodes,historyIndex), ...
        zHistory(normalFaceNodes,historyIndex)];
    arrowDirection = cross( ...
        faceCoordinates(2,:)-faceCoordinates(1,:), ...
        faceCoordinates(3,:)-faceCoordinates(1,:));
    arrowDirection = padding*arrowDirection/norm(arrowDirection);
    set(tractionHandle, ...
        'XData',xHistory(dragNode,historyIndex), ...
        'YData',yHistory(dragNode,historyIndex), ...
        'ZData',zHistory(dragNode,historyIndex), ...
        'UData',arrowDirection(1), ...
        'VData',arrowDirection(2), ...
        'WData',arrowDirection(3));
    title(animationAxes,sprintf('Node %d: step %d/%d', ...
        dragNode,historyIndex-1,completedSteps));
    drawnow;
end
