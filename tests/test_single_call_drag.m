function tests = test_single_call_drag
% Verify one LM call per attempted drag step, including rejected candidates.
tests = functiontests(localfunctions);
end

function testNormalStepsUseOneCallAndAngularIncrement(testCase)
state = runIsolatedDemo('normal',4);
verifyEqual(testCase,state.attemptedSteps,4);
verifyEqual(testCase,state.completedSteps,4);
verifyEqual(testCase,state.totalLmCalls,4);
verifyEqual(testCase,state.totalLmIterations,12);
verifyEqual(testCase,size(state.targets,2),4);
verifyEqual(testCase,size(state.xHistory,2),5);
verifyEqual(testCase,state.failedCorrections,0);
verifyEqual(testCase,state.angleLimitTrials,0);
verifyFalse(testCase,state.angleLimitReached);
verifyEqual(testCase,state.perturbationBisections,0);
verifyFalse(testCase,state.solverInfo.constraintSatisfied);
verifyEqual(testCase,state.maximumResidual,2e-8);

% The mock accepts its target unchanged, exposing each prescribed rotation.
increments = diff([zeros(size(state.targets,1),1),state.targets],1,2);
otherDofs = setdiff(1:size(increments,1),state.dragDofIndices);
verifyEqual(testCase,increments(otherDofs,:),zeros(numel(otherDofs),4));
verifyEqual(testCase,state.dragAngleStep,deg2rad(1));
verifyRotationTargets(testCase,state,ones(1,4));
verifyTiming(testCase,state);
end

function testFirstFailureKeepsInitialStateAndDefinedAverage(testCase)
state = runIsolatedDemo('failure',4);
verifyEqual(testCase,state.attemptedSteps,1);
verifyEqual(testCase,state.completedSteps,0);
verifyEqual(testCase,state.totalLmCalls,1);
verifyEqual(testCase,state.totalLmIterations,3);
verifyEqual(testCase,size(state.targets,2),1);
verifyEqual(testCase,state.failedCorrections,1);
verifyEqual(testCase,state.angleLimitTrials,0);
verifyEqual(testCase,state.xHistory,state.initialCoordinates(:,1));
verifyEqual(testCase,state.yHistory,state.initialCoordinates(:,2));
verifyEqual(testCase,state.zHistory,state.initialCoordinates(:,3));
verifyEqual(testCase,state.maximumResidual,0);
verifyTiming(testCase,state);
end

function testRelativeTargetUsesCorrectedStateEachStep(testCase)
% A correction that holds the initial state must not accumulate target angles.
state = runIsolatedDemo('held-state',3);
verifyEqual(testCase,state.attemptedSteps,3);
verifyEqual(testCase,state.completedSteps,3);
verifyEqual(testCase,state.totalLmCalls,3);
verifyEqual(testCase,size(state.targets,2),3);
verifyGreaterThan(testCase,norm(state.targets(:,1)),0);
verifyEqual(testCase,state.targets,repmat(state.targets(:,1),1,3),'AbsTol',1e-14);
verifyEqual(testCase,state.xHistory,repmat(state.initialCoordinates(:,1),1,4));
verifyEqual(testCase,state.yHistory,repmat(state.initialCoordinates(:,2),1,4));
verifyEqual(testCase,state.zHistory,repmat(state.initialCoordinates(:,3),1,4));
verifyTiming(testCase,state);
end

function testRejectedAttemptPreservesStateAndUsesOneLmCall(testCase)
state = runIsolatedDemo('angle',1);
verifyEqual(testCase,state.attemptedSteps,1);
verifyEqual(testCase,state.completedSteps,0);
verifyEqual(testCase,state.totalLmCalls,1);
verifyEqual(testCase,size(state.targets,2),1);
verifyEqual(testCase,state.failedCorrections,0);
verifyEqual(testCase,state.angleLimitTrials,1);
verifyFalse(testCase,state.angleLimitReached);
verifyEqual(testCase,state.xHistory,state.initialCoordinates(:,1));
verifyEqual(testCase,state.yHistory,state.initialCoordinates(:,2));
verifyEqual(testCase,state.zHistory,state.initialCoordinates(:,3));
verifyEqual(testCase,state.maximumResidual,0);
verifyEqual(testCase,state.maximumAcceptedRotationAngle,deg2rad(170),'AbsTol',1e-14);
verifyTiming(testCase,state);
end

function testOvershootHalvesNextTargetAndKeepsHistoryAcceptedOnly(testCase)
state = runIsolatedDemo('overshoot-recovery',2);
verifyEqual(testCase,state.attemptedSteps,2);
verifyEqual(testCase,state.completedSteps,1);
verifyEqual(testCase,state.totalLmCalls,2);
verifyEqual(testCase,state.totalLmIterations,6);
verifyEqual(testCase,size(state.targets,2),2);
verifyEqual(testCase,state.failedCorrections,0);
verifyEqual(testCase,state.angleLimitTrials,1);
verifyTrue(testCase,state.angleLimitReached);
verifyEqual(testCase,state.perturbationFraction,0.5);
verifyEqual(testCase,state.targets(:,2),0.5*state.targets(:,1),'AbsTol',1e-14);
verifyEqual(testCase,size(state.xHistory,2),2);
verifyEqual(testCase,[state.xHistory(:,1),state.yHistory(:,1),state.zHistory(:,1)], ...
    state.initialCoordinates);
finalDragPosition = [state.xHistory(state.dragNode,end), ...
    state.yHistory(state.dragNode,end),state.zHistory(state.dragNode,end)];
verifyEqual(testCase,finalDragPosition,state.initialCoordinates(state.dragNode,:) ...
    +state.targets(state.dragDofIndices,2).','AbsTol',1e-14);
verifyTiming(testCase,state);
end

function testUnchangedCorrectionStopsBeforeAngleGoal(testCase)
state = runIsolatedDemo('stagnation',4);
verifyEqual(testCase,state.attemptedSteps,1);
verifyEqual(testCase,state.completedSteps,0);
verifyEqual(testCase,state.totalLmCalls,1);
verifyEqual(testCase,size(state.targets,2),1);
verifyGreaterThan(testCase,state.solverInfo.exitflag,0);
verifyEqual(testCase,state.failedCorrections,1);
verifyEqual(testCase,state.angleLimitTrials,0);
verifyFalse(testCase,state.angleLimitReached);
verifyEqual(testCase,state.xHistory,state.initialCoordinates(:,1));
verifyEqual(testCase,state.yHistory,state.initialCoordinates(:,2));
verifyEqual(testCase,state.zHistory,state.initialCoordinates(:,3));
verifyTiming(testCase,state);
end

function testRealLmRunsTwoStepsDespiteRawFeasibilityDiagnostic(testCase)
% A mild increment isolates call semantics from the full default trajectory.
state = runIsolatedDemo('real',2);
verifyEqual(testCase,state.attemptedSteps,2);
verifyEqual(testCase,state.completedSteps,2);
verifyEqual(testCase,state.totalLmCalls,2);
verifyEqual(testCase,state.failedCorrections,0);
verifyEqual(testCase,state.angleLimitTrials,0);
verifyEqual(testCase,size(state.xHistory,2),3);
verifyGreaterThan(testCase,state.solverInfo.exitflag,0);
verifyFalse(testCase,state.solverInfo.constraintSatisfied);
verifyGreaterThan(testCase,state.maximumResidual,1e-30);
verifyTrue(testCase,all(isfinite([state.xHistory(:);state.yHistory(:);state.zHistory(:)])));
verifyTiming(testCase,state);
end

function testNonfiniteCorrectedAngleKeepsLastValidState(testCase)
state = runIsolatedDemo('nonfinite-angle',4);
verifyEqual(testCase,state.attemptedSteps,1);
verifyEqual(testCase,state.completedSteps,0);
verifyEqual(testCase,state.totalLmCalls,1);
verifyEqual(testCase,size(state.targets,2),1);
verifyEqual(testCase,state.failedCorrections,1);
verifyEqual(testCase,state.angleLimitTrials,0);
verifyFalse(testCase,state.angleLimitReached);
verifyEqual(testCase,state.perturbationBisections,0);
verifyEqual(testCase,state.xHistory,state.initialCoordinates(:,1));
verifyEqual(testCase,state.yHistory,state.initialCoordinates(:,2));
verifyEqual(testCase,state.zHistory,state.initialCoordinates(:,3));
verifyEqual(testCase,state.maximumResidual,0);
verifyEqual(testCase,state.maximumAcceptedRotationAngle,0);
verifyTiming(testCase,state);
end

function testStopsNaturallyAtAngleGoal(testCase)
% Run beyond both the former 50-step cap and initial 64-frame allocation.
state = runIsolatedDemo('terminal',66);
verifyEqual(testCase,state.attemptedSteps,66);
verifyEqual(testCase,state.completedSteps,66);
verifyEqual(testCase,state.totalLmCalls,66);
verifyEqual(testCase,size(state.targets,2),66);
verifyEqual(testCase,state.failedCorrections,0);
verifyEqual(testCase,state.angleLimitTrials,0);
verifyTrue(testCase,state.angleLimitReached);
verifyLessThanOrEqual(testCase, ...
    abs(state.maximumAcceptedRotationAngle-pi),state.angleCompletionTolerance);
verifyEqual(testCase,size(state.xHistory,2),67);
verifyTiming(testCase,state);
end

function testInitiallyTerminalStateMakesNoLmCall(testCase)
state = runIsolatedDemo('initial-terminal',0);
verifyEqual(testCase,state.attemptedSteps,0);
verifyEqual(testCase,state.completedSteps,0);
verifyEqual(testCase,state.totalLmCalls,0);
verifyEmpty(testCase,state.targets);
verifyTrue(testCase,state.angleLimitReached);
verifyEqual(testCase,state.xHistory,state.initialCoordinates(:,1));
verifyEqual(testCase,state.yHistory,state.initialCoordinates(:,2));
verifyEqual(testCase,state.zHistory,state.initialCoordinates(:,3));
verifyTiming(testCase,state);
end

function testPredictorBisectsCoordinatePerturbationWithOneLmCall(testCase)
% Bisection scales the Cartesian chord toward the full 1-degree target.
state = runIsolatedDemo('predictor',1);
verifyEqual(testCase,state.attemptedSteps,1);
verifyEqual(testCase,state.completedSteps,1);
verifyEqual(testCase,state.totalLmCalls,1);
verifyEqual(testCase,size(state.targets,2),1);
verifyGreaterThan(testCase,state.perturbationBisections,0);
verifyGreaterThan(testCase,state.perturbationFraction,0);
verifyLessThan(testCase,state.perturbationFraction,1);
verifyLessThanOrEqual(testCase,abs(state.perturbationFraction-0.5),0.005+1e-12);
axisStart = state.initialCoordinates(state.dragRotationEdge(1),:);
axisDirection = state.initialCoordinates(state.dragRotationEdge(2),:)-axisStart;
axisDirection = axisDirection/norm(axisDirection);
initialOffset = state.initialCoordinates(state.dragNode,:)-axisStart;
radial = initialOffset-dot(initialOffset,axisDirection)*axisDirection;
fullPerturbation = (cos(state.dragAngleStep)-1)*radial ...
    +sin(state.dragAngleStep)*cross(radial,axisDirection);
actualPerturbation = state.targets(state.dragDofIndices,1).';
verifyEqual(testCase,actualPerturbation, ...
    state.perturbationFraction*fullPerturbation,'AbsTol',1e-12);
verifyLessThanOrEqual(testCase,state.maximumAcceptedRotationAngle,pi+deg2rad(1e-8));
verifyEqual(testCase,size(state.xHistory,2),2);
verifyTiming(testCase,state);
end

function verifyRotationTargets(testCase,state,fractions)
axisStart = state.initialCoordinates(state.dragRotationEdge(1),:);
axisDirection = state.initialCoordinates(state.dragRotationEdge(2),:)-axisStart;
axisDirection = axisDirection/norm(axisDirection);
positions = [state.initialCoordinates(state.dragNode,:); ...
    state.initialCoordinates(state.dragNode,:)+state.targets(state.dragDofIndices,:).'];
for targetIndex = 1:size(state.targets,2)
    previousOffset = positions(targetIndex,:)-axisStart;
    currentOffset = positions(targetIndex+1,:)-axisStart;
    previousAxial = dot(previousOffset,axisDirection);
    currentAxial = dot(currentOffset,axisDirection);
    previousRadial = previousOffset-previousAxial*axisDirection;
    currentRadial = currentOffset-currentAxial*axisDirection;
    verifyEqual(testCase,currentAxial,previousAxial,'AbsTol',1e-12);
    verifyEqual(testCase,norm(currentRadial),norm(previousRadial),'AbsTol',1e-12);
    panelNormal = cross(previousOffset,axisDirection);
    panelNormal = panelNormal/norm(panelNormal);
    verifyGreaterThan(testCase,dot(currentOffset-previousOffset,panelNormal),0);
    rotation = atan2(norm(cross(previousRadial,currentRadial)), ...
        dot(previousRadial,currentRadial));
    verifyEqual(testCase,rotation,state.dragAngleStep*fractions(targetIndex),'AbsTol',1e-12);
end
end

function verifyTiming(testCase,state)
verifyEqual(testCase,numel(state.stepWallTime),state.attemptedSteps);
verifyTrue(testCase,all(isfinite(state.stepWallTime)));
if state.attemptedSteps == 0
    verifyTrue(testCase,isnan(state.averageStepWallTime));
    return
end
verifyTrue(testCase,isfinite(state.averageStepWallTime));
verifyGreaterThanOrEqual(testCase,state.averageStepWallTime,0);
verifyEqual(testCase,state.averageStepWallTime, ...
    state.solveTime/state.attemptedSteps,'AbsTol',eps);
end

function state = runIsolatedDemo(mode,nSteps)
% Keep mocks, rewritten settings and generated files outside the repository.
codeRoot = fileparts(fileparts(mfilename('fullpath')));
sandboxDirectory = tempname;
mkdir(sandboxDirectory);
originalFolder = pwd;
originalPath = path;
sandboxCleanup = onCleanup(@() restoreSandbox( ...
    originalFolder,originalPath,sandboxDirectory));
files = {'AffineMetricConstraintFactors.m','AffineMetricResidualJacobian.m', ...
    'EvaluateLowRankConstraints.m','MiuraPerturbCorrect.m', ...
    'RotationAngle.m','four_panel_reference.mat'};
for fileIndex = 1:numel(files)
    copyfile(fullfile(codeRoot,files{fileIndex}),sandboxDirectory);
end
scriptText = fileread(fullfile(codeRoot,'demo_miura_drag.m'));
animationStart = strfind(scriptText,'%% Animation');
assert(isscalar(animationStart),'The demo must have one animation section.');
scriptText = scriptText(1:animationStart-1);
scriptText = regexprep(scriptText,'(?m)^(clear|clc|close all);[^\r\n]*(\r?\n|$)','');
if ~strcmp(mode,'terminal') && ~strcmp(mode,'initial-terminal')
    % Bound only the isolated unit-test scenario, never the production loop.
    timerMarker = '    stepTimer = tic;';
    assert(isscalar(strfind(scriptText,timerMarker)),'The step timer was not found.');
    scriptText = strrep(scriptText,timerMarker, ...
        [sprintf('    if attemptedSteps >= %d, break; end\n',nSteps),timerMarker]);
end
if strcmp(mode,'held-state')
    stagnationGuard = 'if ~angleLimitReached && isequal(trialMasterDisplacement,masterDisplacement)';
    assert(contains(scriptText,stagnationGuard),'The stagnation guard was not found.');
    scriptText = strrep(scriptText,stagnationGuard, ...
        ['if false && ',stagnationGuard(4:end)]);
end
if strcmp(mode,'real')
    constraintSetting = '''ConstraintTolerance'',1e-8';
    assert(contains(scriptText,constraintSetting),'The raw feasibility setting was not found.');
    scriptText = strrep(scriptText,constraintSetting,'''ConstraintTolerance'',1e-30');
    incrementSetting = '(?m)^dragAngleStep\s*=\s*deg2rad\(1\)\s*;';
    assert(isscalar(regexp(scriptText,incrementSetting)),'The angular increment was not found.');
    scriptText = regexprep(scriptText,incrementSetting,'dragAngleStep = deg2rad(0.1);');
else
    expectedCalls = nSteps;
    exitFlag = 2;
    if strcmp(mode,'failure') || strcmp(mode,'angle') || ...
            strcmp(mode,'nonfinite-angle') || strcmp(mode,'stagnation')
        expectedCalls = 1;
    end
    if strcmp(mode,'failure')
        exitFlag = -1;
    end
    mockSolver = sprintf([ ...
        'function [beta,info] = MiuraPerturbCorrect(target,varargin)\n' ...
        'logPath = fullfile(fileparts(mfilename(''fullpath'')),''mock_targets.mat'');\n' ...
        'if isfile(logPath), data=load(logPath,''targets''); targets=data.targets;\n' ...
        'else, targets=zeros(numel(target),0); end\n' ...
        'targets(:,end+1)=target; save(logPath,''targets'');\n' ...
        'assert(size(targets,2)<=%d,''Unexpected repeated LM call.'');\n' ...
        'beta=target;\n' ...
        'info=struct(''exitflag'',%d,''iterations'',3,''maxResidual'',2e-8, ...\n' ...
        '    ''constraintSatisfied'',false);\n' ...
        'end\n'],expectedCalls,exitFlag);
    if strcmp(mode,'held-state') || strcmp(mode,'stagnation')
        mockSolver = strrep(mockSolver,'beta=target;','beta=zeros(size(target));');
    end
    writeText(fullfile(sandboxDirectory,'MiuraPerturbCorrect.m'),mockSolver);
    if strcmp(mode,'angle')
        % The predictor stays valid; only the corrected trial crosses the limit.
        mockAngle = [ ...
            'function angles = RotationAngle(nAngles,varargin)' newline ...
            'logPath=fullfile(fileparts(mfilename(''fullpath'')),''mock_targets.mat'');' newline ...
            'if isfile(logPath), angle=deg2rad(-170); else, angle=deg2rad(170); end' newline ...
            'angles=repmat(angle,nAngles,1);' newline 'end' newline];
    elseif strcmp(mode,'overshoot-recovery')
        mockAngle = [ ...
            'function angles = RotationAngle(nAngles,varargin)' newline ...
            'persistent rejectedFirst;' newline ...
            'logPath=fullfile(fileparts(mfilename(''fullpath'')),''mock_targets.mat'');' newline ...
            'calls=0; if isfile(logPath), data=load(logPath,''targets''); calls=size(data.targets,2); end' newline ...
            'if calls>=2, angle=pi;' newline ...
            'elseif calls==1 && isempty(rejectedFirst), angle=deg2rad(-170); rejectedFirst=true;' newline ...
            'else, angle=deg2rad(170); end' newline ...
            'angles=repmat(angle,nAngles,1);' newline 'end' newline];
    elseif strcmp(mode,'nonfinite-angle')
        mockAngle = [ ...
            'function angles = RotationAngle(nAngles,varargin)' newline ...
            'logPath=fullfile(fileparts(mfilename(''fullpath'')),''mock_targets.mat'');' newline ...
            'if isfile(logPath), angle=NaN; else, angle=0; end' newline ...
            'angles=repmat(angle,nAngles,1);' newline 'end' newline];
    elseif strcmp(mode,'initial-terminal')
        mockAngle = ['function angles = RotationAngle(nAngles,varargin)' newline ...
            'angles=pi*ones(nAngles,1);' newline 'end' newline];
    elseif strcmp(mode,'predictor') || strcmp(mode,'terminal')
        mockAngle = [ ...
            'function angles = RotationAngle(nAngles,~,x,y,z)' newline ...
            'persistent initialRadial;' newline ...
            'position=[x(18),y(18),z(18)];' newline ...
            'axisStart=[x(17),y(17),z(17)];' newline ...
            'axisDirection=[x(24),y(24),z(24)]-axisStart;' newline ...
            'axisDirection=axisDirection/norm(axisDirection);' newline ...
            'offset=position-axisStart; radial=offset-dot(offset,axisDirection)*axisDirection;' newline ...
            'if isempty(initialRadial), initialRadial=radial; end' newline ...
            'angleDegree=rad2deg(atan2(norm(cross(initialRadial,radial)),dot(initialRadial,radial)));' newline ...
            'angle=deg2rad(170+20*angleDegree);' newline ...
            'angles=repmat(atan2(sin(angle),cos(angle)),nAngles,1);' newline ...
            'end' newline];
        if strcmp(mode,'terminal')
            mockAngle = strrep(mockAngle,'170+20*angleDegree','114+angleDegree');
        end
    else
        mockAngle = ['function angles = RotationAngle(nAngles,varargin)' newline ...
            'angles=zeros(nAngles,1);' newline 'end' newline];
    end
    writeText(fullfile(sandboxDirectory,'RotationAngle.m'),mockAngle);
end
writeText(fullfile(sandboxDirectory,'demo_miura_drag.m'),scriptText);
addpath(sandboxDirectory,'-begin');
cd(sandboxDirectory);
clear MiuraPerturbCorrect RotationAngle demo_miura_drag
evalc('demo_miura_drag;');

state = struct('attemptedSteps',attemptedSteps,'completedSteps',completedSteps, ...
    'totalLmCalls',totalLmCalls,'totalLmIterations',totalLmIterations, ...
    'failedCorrections',failedCorrections,'angleLimitTrials',angleLimitTrials, ...
    'angleLimitReached',angleLimitReached,'maximumResidual',maximumResidual, ...
    'perturbationBisections',perturbationBisections, ...
    'perturbationFraction',perturbationFraction, ...
    'maximumAcceptedRotationAngle',maximumAcceptedRotationAngle, ...
    'stepWallTime',stepWallTime,'averageStepWallTime',averageStepWallTime, ...
    'solveTime',solveTime,'xHistory',xHistory,'yHistory',yHistory,'zHistory',zHistory, ...
    'initialCoordinates',[x0,y0,z0],'dragDofIndices',dragDofIndices, ...
    'dragNode',dragNode,'dragRotationEdge',dragRotationEdge,'dragAngleStep',dragAngleStep, ...
    'angleCompletionTolerance',angleCompletionTolerance);
if exist('solverInfo','var')
    state.solverInfo = solverInfo;
end
if ~strcmp(mode,'real')
    targetLog = fullfile(sandboxDirectory,'mock_targets.mat');
    if isfile(targetLog)
        targetData = load(targetLog,'targets');
        state.targets = targetData.targets;
    else
        state.targets = zeros(3*nFreeMasters,0);
    end
end
end

function writeText(filePath,text)
fileId = fopen(filePath,'w');
assert(fileId >= 0,'Cannot create an isolated test fixture.');
fileCleanup = onCleanup(@() fclose(fileId));
fprintf(fileId,'%s',text);
end

function restoreSandbox(originalFolder,originalPath,sandboxDirectory)
cd(originalFolder);
path(originalPath);
clear MiuraPerturbCorrect RotationAngle demo_miura_drag
expectedParent = fileparts(fullfile(tempdir,'test_fixture'));
assert(strcmpi(fileparts(sandboxDirectory),expectedParent), ...
    'Refusing to remove a fixture outside the temporary directory.');
if isfolder(sandboxDirectory)
    rmdir(sandboxDirectory,'s');
end
end
