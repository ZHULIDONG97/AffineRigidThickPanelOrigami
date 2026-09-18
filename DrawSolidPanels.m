function h = DrawSolidPanels(x,y,z,panels,varargin)
% Draw semi-transparent solid panels from a cell array of panel faces.

defaultPanelColors = [
    0.16 0.58 0.84
    0.96 0.58 0.22
    0.25 0.66 0.47
    0.58 0.45 0.78
    0.86 0.74 0.24
    ];

parser = inputParser;
parser.addParameter('FaceAlpha',0.42,@(v) isnumeric(v) && isscalar(v) && v >= 0 && v <= 1);
parser.addParameter('EdgeColor',[0.08 0.08 0.08],@(v) isnumeric(v) && numel(v) == 3);
parser.addParameter('LineWidth',1.2,@(v) isnumeric(v) && isscalar(v) && v >= 0);
parser.addParameter('PanelColors',[],@(v) isempty(v) || (isnumeric(v) && size(v,2) == 3));
parser.addParameter('UseLighting',true,@(v) islogical(v) && isscalar(v));
parser.parse(varargin{:});

faceAlpha = parser.Results.FaceAlpha;
edgeColor = parser.Results.EdgeColor;
lineWidth = parser.Results.LineWidth;
panelColors = parser.Results.PanelColors;
useLighting = parser.Results.UseLighting;

x = x(:);
y = y(:);
z = z(:);
nv = numel(x);

if numel(y) ~= nv || numel(z) ~= nv
    error('DrawSolidPanels:CoordinateSizeMismatch', ...
        'x, y, and z must contain the same number of nodes.');
end
if ~iscell(panels)
    error('DrawSolidPanels:InvalidPanels', ...
        'panels must be a cell array.');
end

nPanels = numel(panels);
if isempty(panelColors)
    panelColors = defaultPanelColors;
end

nFacesTotal = 0;
for panelIndex = 1:nPanels
    panelFaces = panels{panelIndex};
    if isempty(panelFaces)
        continue
    elseif iscell(panelFaces)
        nFacesTotal = nFacesTotal + numel(panelFaces);
    elseif isnumeric(panelFaces)
        if isvector(panelFaces)
            nFacesTotal = nFacesTotal + 1;
        else
            nFacesTotal = nFacesTotal + size(panelFaces,1);
        end
    else
        error('DrawSolidPanels:InvalidPanelFaces', ...
            'panels{%d} must be a cell array or numeric face array.', panelIndex);
    end
end

ax = gca;
holdState = ishold(ax);
hold(ax,'on');
h = gobjects(nFacesTotal,1);
faceCounter = 0;

for panelIndex = 1:nPanels
    panelFaces = panels{panelIndex};
    if isempty(panelFaces)
        continue
    end

    colorIndex = mod(panelIndex-1,size(panelColors,1)) + 1;
    faceColor = panelColors(colorIndex,:);

    if iscell(panelFaces)
        for faceIndex = 1:numel(panelFaces)
            face = panelFaces{faceIndex};
            face = face(:).';
            if numel(face) < 3
                error('DrawSolidPanels:InvalidFace', ...
                    'Panel %d face %d must contain at least three nodes.', panelIndex, faceIndex);
            end
            if any(face < 1) || any(face > nv) || any(face ~= round(face))
                error('DrawSolidPanels:InvalidFaceIndex', ...
                    'Panel %d face %d contains an invalid node index.', panelIndex, faceIndex);
            end
            hp = patch(ax, ...
                'XData',x(face), ...
                'YData',y(face), ...
                'ZData',z(face), ...
                'FaceColor',faceColor, ...
                'FaceAlpha',faceAlpha, ...
                'EdgeColor',edgeColor, ...
                'LineWidth',lineWidth, ...
                'FaceLighting','gouraud', ...
                'AmbientStrength',0.45, ...
                'DiffuseStrength',0.72, ...
                'SpecularStrength',0.12, ...
                'SpecularExponent',12);
            faceCounter = faceCounter + 1;
            h(faceCounter,1) = hp;
        end
    elseif isnumeric(panelFaces)
        if isvector(panelFaces)
            faceRows = panelFaces(:).';
        else
            faceRows = panelFaces;
        end
        for faceIndex = 1:size(faceRows,1)
            face = faceRows(faceIndex,:);
            face = face(~isnan(face));
            if numel(face) < 3
                error('DrawSolidPanels:InvalidFace', ...
                    'Panel %d face %d must contain at least three nodes.', panelIndex, faceIndex);
            end
            if any(face < 1) || any(face > nv) || any(face ~= round(face))
                error('DrawSolidPanels:InvalidFaceIndex', ...
                    'Panel %d face %d contains an invalid node index.', panelIndex, faceIndex);
            end
            hp = patch(ax, ...
                'XData',x(face), ...
                'YData',y(face), ...
                'ZData',z(face), ...
                'FaceColor',faceColor, ...
                'FaceAlpha',faceAlpha, ...
                'EdgeColor',edgeColor, ...
                'LineWidth',lineWidth, ...
                'FaceLighting','gouraud', ...
                'AmbientStrength',0.45, ...
                'DiffuseStrength',0.72, ...
                'SpecularStrength',0.12, ...
                'SpecularExponent',12);
            faceCounter = faceCounter + 1;
            h(faceCounter,1) = hp;
        end
    else
        error('DrawSolidPanels:InvalidPanelFaces', ...
            'panels{%d} must be a cell array or numeric face array.', panelIndex);
    end
end

h = h(1:faceCounter);
if useLighting && ~isempty(h)
    delete(findall(ax,'Type','light','Tag','DrawSolidPanelsLight'));
    light(ax,'Position',[0.45 -0.60 0.75], ...
        'Style','infinite', ...
        'Color',[1 1 1], ...
        'Tag','DrawSolidPanelsLight');
    lighting(ax,'gouraud');
end
axis(ax,'equal');
if ~holdState
    hold(ax,'off');
end
end
