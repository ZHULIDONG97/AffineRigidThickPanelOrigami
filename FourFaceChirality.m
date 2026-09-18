close all
clear
clc

ra=pi/180;
Ra=1; %radius
constraintRankTolerance = 1e-10;

Rf=[60*ra 60*ra 60*ra 60*ra 60*ra 60*ra];

n=length(Rf);
theta=zeros(n,1);
for i=1:n
    if i==1
        theta(i)=0;
    else
        theta_temp=0;
        for j=1:i-1
            theta_temp=theta_temp+Rf(j);
        end
        theta(i)=theta_temp;
    end
end

% 定义顶点坐标（三维坐标，z 坐标为 0）
Xcorner=zeros(n,1);
Ycorner=zeros(n,1);
Zcorner=zeros(n,1);
for i=1:n
Xcorner(i)=Ra*cos(theta(i));
Ycorner(i)=Ra*sin(theta(i));
Zcorner(i)=0;
end
XYZcorner=[Xcorner Ycorner Zcorner];

% figure
% scatter(Xcorner,Ycorner);
% axis off; axis equal; axis tight;

vertices = [Xcorner(1) Ycorner(1) 0; %1
            Xcorner(2) Ycorner(2) 0; %2
            Xcorner(2) Ycorner(2) 0; %3
            Xcorner(3) Ycorner(3) 0;%4
            Xcorner(4) Ycorner(4) 0;%5
            Xcorner(4) Ycorner(4) 0;%6
            Xcorner(5) Ycorner(5) 0;%7
            Xcorner(6) Ycorner(6) 0;%8
            Xcorner(6) Ycorner(6) 0;%9
            Xcorner(1) Ycorner(1) 0.1;%10
            Xcorner(1) Ycorner(1) 0.1;%11
            Xcorner(2) Ycorner(2) 0.1;%12
            Xcorner(3) Ycorner(3) 0.1;%13
            Xcorner(5) Ycorner(5) 0.1;%14
            Xcorner(6) Ycorner(6) 0.1;%15
            Xcorner(3) Ycorner(3) 0.2;%16
            Xcorner(4) Ycorner(4) 0.2;%17
            Xcorner(5) Ycorner(5) 0.2;%18
            0 0 0;%19
            0 0 0.1;%20
            0 0 0.1;%21
            0 0 0;%22
            0 0 0;%23
            0 0 0.2;%24
            ];

edges1=[1 2; 2 19; 19 1;
        11 12; 12 21; 21 11;
        1 11; 2 12; 19 21
        ];

edges2=[1 19; 19 9;9 1;
        10 20; 20 15; 15 10;
        1 10; 19 20; 9 15
        ];

edges3=[3 4; 4 5; 5 23; 23 3;
        12 13; 16 17; 17 24; 21 12;
        13 21; 16 24;
        3 12; 4 13; 13 16; 5 17; 23 21; 21 24
        ];

edges4=[6 7; 7 8; 8 22; 22 6;
        17 18; 14 15; 15 20; 24 17;
        14 20; 18 24;
        6 17; 7 14; 14 18; 8 15; 22 20; 20 24
       ];

edges = [edges1;edges2;edges3;edges4;];

panels = cell(4,1);

panels{1} = {[1 2 19],[11 12 21],[1 11 12 2],[2 12 21 19],[1 11 21 19]};
panels{2} = {[1 19 9],[10 20 15],[1 10 15 9],[1 10 20 19],[9 15 20 19]};
panels{3} = {[5 23 3 4],[12 13 21],[16 17 24],[5 23 24 17],[3 12 21 23],[3 12 13 4],[4 16 17 5],[13 16 24 21]};
panels{4} = {[6 7 8 22],[17 18 24],[14 15 20],[6 7 18 17],[7 8 15 14],[8 22 20 15],[6 22 24 17],[14 20 24 18]};

% define number and sequence
nv=size(vertices,1);
ne=size(edges,1);


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
% face3
p = unique(edges3);
m = numel(p);
B_e = [ones(m,1) x0(p) y0(p) z0(p)];
E_e = eye(m) - B_e * pinv(B_e,1e-5);
E(p, p) = E(p, p) + E_e;
% face4
p = unique(edges4);
m = numel(p);
B_e = [ones(m,1) x0(p) y0(p) z0(p)];
E_e = eye(m) - B_e * pinv(B_e,1e-5);
E(p, p) = E(p, p) + E_e;

% Use the validated global master set to parameterize the affine space.
masterNodeIndices = [5 16 17 24 6 18 1 19];
[T, r] = TransformMain(masterNodeIndices, nv);

Em = E(:,masterNodeIndices);
Es = E(:,setdiff(1:nv,masterNodeIndices));

svd(Es)

B = [eye(r);-pinv(Es,1e-5)*Em];
mappedB = T \ B;

%% boundary condition
fix = [16 17 24];
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
Cchi=[19 1 2 11; 9 1 19 10; 24 17 16 5; 18 17 24 6];

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
Cfold = [2 19 1 9 -1];

% rng(1,'twister');
perturbScale = 1;
chiMargin = 0;
foldAngleMargin = 0;
equalityPenaltySchedule = 1e6;

betaTarget = perturbScale * randn(3*r,1);
beta0 = betaTarget;

[cineq0, ceq0] = ConstraintsGen(beta0,x0,y0,z0,B, ...
    rigidityDirections,rigidityReferences,rigidityOffsets, ...
    Cchi,Cfold,r,T,chiMargin,foldAngleMargin);
fprintf('Initial: ||ceq|| = %.3e, max(cineq) = %.3e\n', norm(ceq0), max(cineq0));

options = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','iter', ...
    'MaxIterations',10000, ...
    'MaxFunctionEvaluations',5e5, ...
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

%% Plotting

% panelColors = [
%     0.16 0.58 0.84
%     0.96 0.58 0.22
%     0.25 0.66 0.47
%     0.58 0.45 0.78
%     ];
panelColors = [
    0.00 0.45 0.85   % vivid blue
    1.00 0.45 0.00   % vivid orange
    0.00 0.70 0.35   % vivid green
    0.75 0.25 0.95   % vivid purple
    ];
panelAlpha = 0.7;
panelEdgeWidth = 2.0;

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

figure
hold on
DrawSolidPanels(x,y,z,panels, ...
    'PanelColors',panelColors, ...
    'FaceAlpha',0.75, ...
    'LineWidth',panelEdgeWidth);
scatter3(x(masterNodeIndices),y(masterNodeIndices),z(masterNodeIndices), ...
    70,'r','filled')
% panelOtherNodes = setdiff(1:nv,mainPointIndices);
% scatter3(x(panelOtherNodes),y(panelOtherNodes),z(panelOtherNodes), ...
%     70,'k','filled')
axis off; axis equal; axis tight;
view([5,30])
