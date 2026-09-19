function tests = test_panel_faces
tests = functiontests(localfunctions);
end

function testMixedFaceRepresentations(testCase)
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(fileparts(fileparts(mfilename('fullpath'))));
fig = figure('Visible','off');
figureCleanup = onCleanup(@() close(fig));
ax = axes(fig);
x = [0;1;1;0];
y = [0;0;1;1];
z = zeros(4,1);
panels = {{1:4},[1 2 3 NaN;1 3 4 NaN],(1:4).',[]};
colors = [1 0 0;0 1 0];
h = DrawSolidPanels(x,y,z,panels,'PanelColors',colors,'FaceAlpha',0.6);
faces = {1:4,[1 2 3],[1 3 4],1:4};
verifySize(testCase,h,[4 1]);
for i = 1:numel(h)
    verifyEqual(testCase,h(i).XData(:),x(faces{i}));
    verifyEqual(testCase,h(i).YData(:),y(faces{i}));
    verifyEqual(testCase,h(i).ZData(:),z(faces{i}));
    verifyEqual(testCase,h(i).FaceAlpha,0.6);
end
verifyEqual(testCase,vertcat(h.FaceColor),colors([1 2 2 1],:));
verifyFalse(testCase,ishold(ax));
hold(ax,'on');
DrawSolidPanels(x,y,z,panels);
verifyTrue(testCase,ishold(ax));
verifyNumElements(testCase,findall(ax,'Type','light','Tag','DrawSolidPanelsLight'),1);
verifyError(testCase,@() DrawSolidPanels(x,y,z,{{[1 2 NaN]}}), ...
    'DrawSolidPanels:InvalidFaceIndex');
verifyError(testCase,@() DrawSolidPanels(x,y,z,{[1 2 NaN]}), ...
    'DrawSolidPanels:InvalidFace');
verifyError(testCase,@() DrawSolidPanels(x,y,z,{[1 2 5]}), ...
    'DrawSolidPanels:InvalidFaceIndex');
verifyError(testCase,@() DrawSolidPanels(x,y,z,{struct()}), ...
    'DrawSolidPanels:InvalidPanelFaces');
end
