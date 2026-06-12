function fig = plot_sloreta_signed_cortex(sourceXYZ, sourceValues, varargin)
%PLOT_SLORETA_SIGNED_CORTEX Plot signed sLORETA contrast values on cortex.
%
%   fig = plot_sloreta_signed_cortex(sourceXYZ, sourceValues)
%
%   Positive values are plotted with warm colors and mean:
%
%       condition 1 > condition 2
%
%   Negative values are plotted with cool colors and mean:
%
%       condition 2 > condition 1
%
%   Name-value options
%   ------------------
%   'View'               'six', 'left', 'top', 'right', 'front', 'back',
%                        'bottom', or [az el]. Default: 'six'.
%   'Threshold'          Values with abs(value) below this are gray.
%                        Default: 0.
%   'ColorLimit'         Symmetric color limit. Default: max abs value.
%   'ProjectionRadius'   Source-to-cortex projection radius in mm.
%                        Default: 35.
%   'ProjectionMode'     '3d' or 'yz'. Default: '3d'.
%   'BrainTemplate'      'brainnet', 'brainnet_smoothed', or 'brainstorm'.
%                        Default: 'brainnet'.
%   'SurfaceFile'        Optional custom .mat or BrainNet .nv surface file.
%   'OutputPng'          Optional PNG path to export.
%   'ShowTitle'          true/false. Default: false.
%   'TitleText'          Optional title text.

parser = inputParser;
parser.FunctionName = mfilename;
addRequired(parser, 'sourceXYZ', @(x) isnumeric(x) && size(x, 2) == 3);
addRequired(parser, 'sourceValues', @(x) isnumeric(x) && isvector(x));
addParameter(parser, 'View', 'six', @(x) (ischar(x) || isstring(x)) || ...
    (isnumeric(x) && numel(x) == 2));
addParameter(parser, 'Threshold', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(parser, 'ColorLimit', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && x > 0));
addParameter(parser, 'ProjectionRadius', 35, @(x) isnumeric(x) && ...
    isscalar(x) && isfinite(x) && x > 0);
addParameter(parser, 'ProjectionMode', '3d', @(x) ischar(x) || isstring(x));
addParameter(parser, 'BrainTemplate', 'brainnet', @(x) ischar(x) || isstring(x));
addParameter(parser, 'SurfaceFile', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputPng', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ShowTitle', false, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'TitleText', '', @(x) ischar(x) || isstring(x));
parse(parser, sourceXYZ, sourceValues, varargin{:});

sourceValues = sourceValues(:);
if size(sourceXYZ, 1) ~= numel(sourceValues)
    error('%s:SizeMismatch', mfilename, ...
        'sourceXYZ rows must match sourceValues length.');
end

[vertices, faces, surfaceKind] = load_cortex_surface(char(parser.Results.BrainTemplate), ...
    char(parser.Results.SurfaceFile));

vertexValues = project_signed_sources_to_vertices(vertices, sourceXYZ, sourceValues, ...
    parser.Results.ProjectionRadius, char(parser.Results.ProjectionMode));
vertexValues = smooth_vertex_values(faces, vertexValues, 3);

threshold = parser.Results.Threshold;
colorLimit = parser.Results.ColorLimit;
if isempty(colorLimit)
    colorLimit = max(abs(vertexValues));
end
if isempty(colorLimit) || colorLimit <= 0
    colorLimit = 1;
end

colors = signed_cortex_colors(vertexValues, threshold, colorLimit);

viewSpec = parser.Results.View;
fig = figure('Color', 'k', 'Name', 'Signed sLORETA cortex map');
set(fig, 'Renderer', 'opengl', 'InvertHardcopy', 'off');

if ischar(viewSpec) || isstring(viewSpec)
    isSixView = strcmpi(char(viewSpec), 'six');
else
    isSixView = false;
end

if isSixView
    set(fig, 'Position', [80 80 1500 900]);
    viewNames = {'left', 'top', 'right', 'front', 'back'};
    axesPositions = [0.05 0.54 0.28 0.36; ...
        0.36 0.54 0.28 0.36; ...
        0.67 0.54 0.28 0.36; ...
        0.20 0.10 0.28 0.36; ...
        0.52 0.10 0.28 0.36];
    for viewIdx = 1:numel(viewNames)
        ax = axes('Parent', fig, 'Position', axesPositions(viewIdx, :));
        render_cortex_axes(ax, vertices, faces, colors, viewNames{viewIdx}, 1.02, surfaceKind);
        title(ax, viewNames{viewIdx}, 'Color', 'w', 'FontWeight', 'normal');
    end
    if parser.Results.ShowTitle
        titleText = char(parser.Results.TitleText);
        annotation(fig, 'textbox', [0.02 0.955 0.96 0.04], ...
            'String', titleText, 'Color', 'w', 'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
    end
else
    set(fig, 'Position', [100 100 1100 820]);
    ax = axes('Parent', fig, 'Color', 'k');
    set(ax, 'Position', [0 0 1 1]);
    render_cortex_axes(ax, vertices, faces, colors, viewSpec, 1.12, surfaceKind);
    if parser.Results.ShowTitle
        title(ax, char(parser.Results.TitleText), 'Color', 'w', 'FontWeight', 'normal');
    end
end

outputPng = char(parser.Results.OutputPng);
if ~isempty(outputPng)
    save_figure_png(fig, outputPng, 250);
end

end

function [vertices, faces, surfaceKind] = load_cortex_surface(brainTemplate, surfaceFile)
rootDir = fileparts(mfilename('fullpath'));
if isempty(surfaceFile)
    surfaceFile = builtin_surface_file(rootDir, brainTemplate);
end
surfaceKind = surface_kind_from_template(brainTemplate, surfaceFile);
if ~exist(surfaceFile, 'file')
    error('%s:MissingSurface', mfilename, ...
        'Could not find cortex surface file: %s', surfaceFile);
end

function surfaceKind = surface_kind_from_template(brainTemplate, surfaceFile)
[~, ~, ext] = fileparts(surfaceFile);
if strcmpi(ext, '.nv') || startsWith(lower(char(brainTemplate)), 'brainnet')
    surfaceKind = 'brainnet';
else
    surfaceKind = 'brainstorm';
end
end

[~, ~, ext] = fileparts(surfaceFile);
switch lower(ext)
    case '.mat'
        cortex = load(surfaceFile, 'Vertices', 'Faces');
        if isfield(cortex, 'Vertices') && isfield(cortex, 'Faces')
            vertices = cortex.Vertices;
            faces = cortex.Faces;
        else
            loaded = load(surfaceFile, 'mesh');
            vertices = loaded.mesh.pos;
            faces = loaded.mesh.tri;
        end
    case '.nv'
        [vertices, faces] = read_brainnet_nv(surfaceFile);
    otherwise
        error('%s:BadSurfaceFile', mfilename, ...
            'SurfaceFile must be a .mat or BrainNet .nv file.');
end

if max(abs(vertices(:))) < 1
    vertices = vertices * 1000;
end
end

function surfaceFile = builtin_surface_file(rootDir, brainTemplate)
templateDir = fullfile(rootDir, 'templates', 'cortex');
switch lower(char(brainTemplate))
    case {'brainnet', 'brainnet_icbm152', 'icbm152'}
        surfaceFile = fullfile(templateDir, 'brainnet_icbm152.nv');
    case {'brainnet_smoothed', 'brainnetsmoothed', 'icbm152_smoothed'}
        surfaceFile = fullfile(templateDir, 'brainnet_icbm152_smoothed.nv');
    case {'brainstorm', 'brainstorm_low', 'brainstorm_pial'}
        surfaceFile = fullfile(templateDir, 'brainstorm_icbm152_cortex_pial_low.mat');
    otherwise
        error('%s:BadBrainTemplate', mfilename, ...
            'BrainTemplate must be brainnet, brainnet_smoothed, or brainstorm.');
end
end

function [vertices, faces] = read_brainnet_nv(surfaceFile)
fid = fopen(surfaceFile, 'r');
if fid < 0
    error('%s:CannotOpenSurface', mfilename, ...
        'Could not open BrainNet surface file: %s', surfaceFile);
end
cleaner = onCleanup(@() fclose(fid));
nVertices = fscanf(fid, '%d', 1);
vertices = fscanf(fid, '%f', [3 nVertices]).';
nFaces = fscanf(fid, '%d', 1);
faces = fscanf(fid, '%d', [3 nFaces]).';
if min(faces(:)) == 0
    faces = faces + 1;
end
end

function vertexValues = project_signed_sources_to_vertices(vertices, sourceXYZ, sourceValues, radius, mode)
activeIdx = find(sourceValues ~= 0);
if isempty(activeIdx)
    vertexValues = zeros(size(vertices, 1), 1);
    return;
end

activeXYZ = sourceXYZ(activeIdx, :);
activeValues = sourceValues(activeIdx);
nVertices = size(vertices, 1);
vertexValues = zeros(nVertices, 1);
weightSums = zeros(nVertices, 1);
chunkSize = 1500;
sigma2 = 2 * radius ^ 2;

switch lower(char(mode))
    case '3d'
        vertexCoords = vertices;
        sourceCoords = activeXYZ;
    case 'yz'
        vertexCoords = vertices(:, 2:3);
        sourceCoords = activeXYZ(:, 2:3);
    otherwise
        error('%s:BadProjectionMode', mfilename, ...
            'ProjectionMode must be ''3d'' or ''yz''.');
end

sourceNorm = sum(sourceCoords .^ 2, 2).';
for firstIdx = 1:chunkSize:nVertices
    lastIdx = min(firstIdx + chunkSize - 1, nVertices);
    v = vertexCoords(firstIdx:lastIdx, :);
    dist2 = sum(v .^ 2, 2) + sourceNorm - 2 * (v * sourceCoords.');
    weights = exp(-dist2 / sigma2);
    vertexValues(firstIdx:lastIdx) = weights * activeValues;
    weightSums(firstIdx:lastIdx) = sum(weights, 2);
end

valid = weightSums > 0;
vertexValues(valid) = vertexValues(valid) ./ weightSums(valid);
end

function values = smooth_vertex_values(faces, values, nIter)
nVertices = numel(values);
edges = [faces(:, [1 2]); faces(:, [2 3]); faces(:, [3 1])];
edges = [edges; edges(:, [2 1])];
for iter = 1:nIter
    accum = accumarray(edges(:, 1), values(edges(:, 2)), [nVertices 1], @sum, 0);
    counts = accumarray(edges(:, 1), 1, [nVertices 1], @sum, 0);
    values = 0.55 * values + 0.45 * (accum ./ max(counts, 1));
end
end

function colors = signed_cortex_colors(values, threshold, colorLimit)
baseGray = repmat([0.58 0.58 0.56], numel(values), 1);
colors = baseGray;

positive = values > threshold;
negative = values < -threshold;

if any(positive)
    scaled = min(values(positive) ./ colorLimit, 1);
    warm = hot(256);
    idx = max(1, min(256, round(1 + scaled * 255)));
    colors(positive, :) = warm(idx, :);
end

if any(negative)
    scaled = min(abs(values(negative)) ./ colorLimit, 1);
    coolMap = cool(256);
    idx = max(1, min(256, round(1 + scaled * 255)));
    colors(negative, :) = coolMap(idx, :);
end
end

function render_cortex_axes(ax, vertices, faces, colors, viewSpec, zoomFactor, surfaceKind)
set(ax, 'Color', 'k');
patch('Parent', ax, ...
    'Vertices', vertices, ...
    'Faces', faces, ...
    'FaceVertexCData', colors, ...
    'FaceColor', 'interp', ...
    'EdgeColor', 'none', ...
    'SpecularStrength', 0.10, ...
    'DiffuseStrength', 0.86, ...
    'AmbientStrength', 0.34);
axis(ax, 'equal');
axis(ax, 'off');
view(ax, view_to_az_el(viewSpec, surfaceKind));
camlight(ax, 'headlight');
camlight(ax, -80, 20);
camlight(ax, 100, 25);
lighting(ax, 'gouraud');
material(ax, 'dull');
camzoom(ax, zoomFactor);
end

function azel = view_to_az_el(viewSpec, surfaceKind)
if isnumeric(viewSpec)
    azel = viewSpec(:).';
    return;
end

if strcmpi(surfaceKind, 'brainnet')
    switch lower(char(viewSpec))
        case 'right'
            azel = [-90 0];
        case 'left'
            azel = [90 0];
        case 'top'
            azel = [0 90];
        case 'front'
            azel = [180 0];
        case 'back'
            azel = [0 0];
        case 'bottom'
            azel = [0 -90];
        otherwise
            error('%s:BadView', mfilename, ...
                'View must be six, left, top, right, front, back, bottom, or [az el].');
    end
else
    switch lower(char(viewSpec))
        case 'right'
            azel = [0 0];
        case 'left'
            azel = [180 0];
        case 'top'
            azel = [0 90];
        case 'front'
            azel = [90 0];
        case 'back'
            azel = [-90 0];
        case 'bottom'
            azel = [0 -90];
        otherwise
            error('%s:BadView', mfilename, ...
                'View must be six, left, top, right, front, back, bottom, or [az el].');
    end
end
end

function save_figure_png(fig, outputPng, resolution)
if exist('exportgraphics', 'file') == 2
    exportgraphics(fig, outputPng, 'Resolution', resolution, ...
        'BackgroundColor', 'black');
else
    set(fig, 'InvertHardcopy', 'off');
    print(fig, outputPng, '-dpng', sprintf('-r%d', resolution));
end
end
