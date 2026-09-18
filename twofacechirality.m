close all
clear
clc

constraintRankTolerance = 1e-10;

vertices = [
    0 0 0;
    1 0 0;
    1 1 0;
    0 1 0;
    0 0 0.2;
    1 0 0.2;
    1 1 0.2;
    0 1 0.2;
    1 0 0.2;
    1 1 0.2;
    2 0 0;
    2 1 0;
    2 0 0.2
    2 1 0.2;
    ];


edges1=[1 2; 2 3; 3 4; 4 1;
       5 6; 6 7; 7 8; 8 5;
       1 5; 2 6; 3 7; 4 8
       ];

edges2=[2 11; 11 12; 12 3; 3 2;
        9 13; 13 14; 14 10; 10 9;
        2 9; 11 13; 12 14; 3 10
        ];

edges = [edges1;edges2];

panels = cell(2,1);

panels{1} = {[1 2 3 4],[5 6 7 8],[1 2 6 5],[4 3 7 8],[2 3 7 6],[1 4 8 5]};
panels{2} = {[2 11 12 3],[9 13 14 10],[2 11 13 9],[3 12 14 10],[11 12 14 13],[2 3 10 9]};

% define number and sequence
nv=size(vertices,1);
ne=size(edges,1);

nf = 2;

CON = edges;
C=zeros(ne,nv);
for i = 1:ne      % connectivity matrix
    C( i, CON(i,1) ) = 1;
    C( i, CON(i,2) ) = -1;
end

x0=vertices(:,1);
y0=vertices(:,2);
z0=vertices(:,3);

u = C*x0;
v = C*y0;
w = C*z0;
l0 = (u.^2 + v.^2 + w.^2) .^0.5;

%% build affine projection matrix
E = zeros(nv,nv);
%face1
p = unique(edges1);
m = numel(p);
B_e = [ones(m,1) x0(p) y0(p) z0(p)];
E_e = eye(m) - B_e * pinv(B_e,1e-5);
E(p, p) = E(p, p) + E_e;
% face2
p = unique(edges2);
m = numel(p);
B_e = [ones(m,1) x0(p) y0(p) z0(p)];
E_e = eye(m) - B_e * pinv(B_e,1e-5);
E(p, p) = E(p, p) + E_e;

mainPointIndices = [1 2 3 6 9 11];
% r = numel(mainPointIndices);
[T, r] = TransformMain(mainPointIndices, nv);

Em = E(:,mainPointIndices);
Es = E(:,setdiff(1:nv,mainPointIndices));

svd(Es)

B = [eye(r);-pinv(Es,1e-5)*Em];
mappedB = T \ B;

%% boundary condition
fix = [1 2 3];
% fix = [];
% Keep fixed supports out of E so the six main-node coordinates remain.
% They are imposed below as hard linear constraints on beta.
if isempty(fix)
    AeqFix = zeros(0,3*r);
    beqFix = zeros(0,1);
else

fixedMap = mappedB(fix,:);
nfix = numel(fix);
AeqFix = zeros(3*nfix,3*r);
AeqFix(1:nfix,1:r) = fixedMap;
AeqFix(nfix+1:2*nfix,r+1:2*r) = fixedMap;
AeqFix(2*nfix+1:3*nfix,2*r+1:3*r) = fixedMap;
beqFix = zeros(3*nfix,1);
end

%% folding path

XYZ0=[x0; y0; z0];

% Orientation rows keep each thick panel from flipping into its mirror image.
Cchi=[1 2 3 6; 3 2 11 9];

% Preset six metric components per solid panel and retain an independent set.
panelNodes = cellfun(@(panelFaces) unique([panelFaces{:}],'stable'), ...
    panels,'UniformOutput',false);
[rigidityDirections,rigidityReferences,rigidityOffsets, ...
    metricConstraints,candidateMetricConstraints] = ...
    AffineMetricConstraintFactors(mappedB,x0,y0,z0,panelNodes,3, ...
    constraintRankTolerance);
fprintf('Thick-panel metric constraints: %d preset, %d independent.\n', ...
    size(candidateMetricConstraints,1),size(metricConstraints,1));

% Valley fold branch on crease edge 2-3.  The sign is defined by
% the four points [1 2 3 12].
Cfold = [1 2 3 12 -1];

% rng(1,'twister');
perturbScale = 1;
chiMargin = 1e-8;
foldAngleMargin = 1e-4;
equalityPenaltySchedule = 1e6;

betaTarget = perturbScale * randn(3*r,1);
beta0 = betaTarget;

[cineq0, ceq0] = ConstraintsGen(beta0,x0,y0,z0,B, ...
    rigidityDirections,rigidityReferences,rigidityOffsets, ...
    Cchi,Cfold,r,T,chiMargin,foldAngleMargin);
fprintf('Initial: ||ceq|| = %.3e, max(cineq) = %.3e\n', norm(ceq0), max(cineq0));

options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','iter', ...
    'MaxIterations',1000, ...
    'MaxFunctionEvaluations',5e4, ...
    'ConstraintTolerance',1e-10, ...
    'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-12, ...
    'FiniteDifferenceType','central');

nonlcon = @(beta) inequalityOnlyConstraint(beta,x0,y0,z0,B, ...
    rigidityDirections,rigidityReferences,rigidityOffsets,Cchi,Cfold, ...
    r,T,chiMargin,foldAngleMargin);
objective = @(beta) equalityPenaltyObjective(beta,betaTarget,equalityPenaltySchedule, ...
        x0,y0,z0,B,rigidityDirections,rigidityReferences, ...
        rigidityOffsets,Cchi,Cfold,r,T,chiMargin,foldAngleMargin);

tic

    [betafs, objectiveValue, exitflag, output] = fmincon(objective,beta0, ...
        [],[],AeqFix,beqFix,[],[],nonlcon,options);
    
CPU=toc;

[cineqFinal, ceqFinal, constraintState] = ConstraintsGen( ...
    betafs,x0,y0,z0,B,rigidityDirections,rigidityReferences, ...
    rigidityOffsets,Cchi,Cfold,r,T,chiMargin,foldAngleMargin);

constraintResidual = [ceqFinal; max(0,cineqFinal)];
disp('Constraint residual [ceq; active inequality violation]:');
disp(constraintResidual);

fprintf('Final: ||ceq|| = %.3e, max(cineq) = %.3e\n', norm(ceqFinal), max(cineqFinal));
fprintf(' fold angle = %.3e rad\n',constraintState.foldAngle);
fprintf('Computational time (s) =');
disp(CPU);

x = constraintState.x;
y = constraintState.y;
z = constraintState.z;
    
u = C*x;
v = C*y;
w = C*z;
l = (u.^2 + v.^2 + w.^2) .^0.5;

lengthError = l-l0;
lengthErrorNorm = norm(lengthError);
% disp('Edge length error:');
% disp(lengthError);
fprintf('Edge length error norm = %.3e\n', lengthErrorNorm);

%% 
close all

panelColors = [
    0.16 0.58 0.84
    0.96 0.58 0.22
    ];
panelAlpha = 0.7;
panelEdgeWidth = 1.5;

figure
hold on
DrawSolidPanels(x0,y0,z0,panels, ...
    'PanelColors',panelColors, ...
    'FaceAlpha',panelAlpha, ...
    'LineWidth',panelEdgeWidth);

for k = 1:size(edges,1)
    idx_start = edges(k, 1);
    idx_end = edges(k, 2);
    plot3([x0(idx_start), x0(idx_end)], [y0(idx_start), y0(idx_end)],...
        [z0(idx_start), z0(idx_end)], 'k-', 'LineWidth', 0.8);
        
    text((x0(idx_start) + x0(idx_end))/2, (y0(idx_start) + y0(idx_end))/2,...
        (z0(idx_start) + z0(idx_end))/2, num2str(k), ...
        'FontSize', 10, ...
        'Color', 'r', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end
for j=1:size(vertices,1)   % node sequence
    text(vertices(j,1), vertices(j,2),vertices(j,3), [' ' num2str(j)], 'FontSize', 12, 'Color', 'b');
end
axis off; axis equal; axis tight;
view([5,30])

figure
hold on
DrawSolidPanels(x0,y0,z0,panels, ...
    'PanelColors',panelColors, ...
    'FaceAlpha',panelAlpha, ...
    'LineWidth',panelEdgeWidth);
scatter3(x0(fix),y0(fix),z0(fix),48,'k','filled')

axis off; axis equal; axis tight;
view([5,30])

figure
hold on
for k = 1:size(edges,1)
    idx_start = edges(k, 1);
    idx_end = edges(k, 2);
    plot3([x0(idx_start), x0(idx_end)], [y0(idx_start), y0(idx_end)],...
        [z0(idx_start), z0(idx_end)], '--','Color',[0.5 0.5 0.5],'LineWidth',1.0);
end
DrawSolidPanels(x,y,z,panels, ...
    'PanelColors',panelColors, ...
    'FaceAlpha',panelAlpha, ...
    'LineWidth',panelEdgeWidth);
axis off; axis equal; axis tight;
view([5,30])

figure
hold on
DrawSolidPanels(x,y,z,panels, ...
    'PanelColors',panelColors, ...
    'FaceAlpha',panelAlpha, ...
    'LineWidth',panelEdgeWidth);
% scatter3(x0(fix),y0(fix),z0(fix),48,'k','filled')
axis off; axis equal; axis tight;
view([5,30])





