function fig = plot_sloreta_cortex(sourceFile, varargin)
%PLOT_SLORETA_CORTEX Plot sLORETA values on a folded cortical surface.
%
%   fig = plot_sloreta_cortex(sourceFile)
%   fig = plot_sloreta_cortex(sourceFile, 'TimeIndex', 317)
%
%   This makes a presentation-style cortical plot using a local copy of the
%   BrainNet ICBM152 pial cortex. Source-grid values are interpolated to
%   the nearest cortical vertices.
%
%   Name-value options
%   ------------------
%   'TimeIndex'          Time sample to plot. Default: peak mean source value.
%   'ThresholdPercentile' Plot warm colors above this percentile. Default: 88.
%   'ColorPercentile'    Upper color-limit percentile. Default: 99.
%   'View'               'left', 'top', 'right', 'front', 'back', 'bottom',
%                        'six', or [az el]. The 'six' option creates a
%                        subplot(2,3) montage with left/top/right on the top
%                        row and front/back/bottom on the bottom row.
%                        Default: 'right'.
%   'Region'             'all' or 'acc'. Default: 'all'.
%   'ACCBox'             ACC bounding box [xmin xmax ymin ymax zmin zmax].
%                        Default: [-15 15 0 50 0 55].
%   'MaxSources'         Keep only the strongest N sources after masking.
%                        Inf keeps all masked sources. Default: Inf.
%   'ProjectionRadius'   Radius in mm for projecting active sources to the
%                        visible cortex. Default: 35.
%   'ProjectionMode'     '3d' or 'yz'. Use 'yz' to display medial ACC on a
%                        lateral brain silhouette. Default: '3d'.
%   'BrainTemplate'      'brainnet', 'brainnet_smoothed', or 'brainstorm'.
%                        Default: 'brainnet'.
%   'SurfaceFile'        Optional custom .mat or BrainNet .nv surface file.
%   'ShowTitle'          true/false. Default: false.
%   'OutputPng'          Optional PNG path to export.

parser = inputParser;
parser.FunctionName = mfilename;
addRequired(parser, 'sourceFile', @(x) ischar(x) || isstring(x));
addParameter(parser, 'TimeIndex', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1));
addParameter(parser, 'ThresholdPercentile', 88, @(x) isnumeric(x) && ...
    isscalar(x) && isfinite(x) && x >= 0 && x <= 100);
addParameter(parser, 'ColorPercentile', 99, @(x) isnumeric(x) && ...
    isscalar(x) && isfinite(x) && x >= 0 && x <= 100);
addParameter(parser, 'View', 'right', @(x) (ischar(x) || isstring(x)) || ...
    (isnumeric(x) && numel(x) == 2));
addParameter(parser, 'Region', 'all', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ACCBox', [-15 15 0 50 0 55], @(x) isnumeric(x) && numel(x) == 6);
addParameter(parser, 'MaxSources', Inf, @(x) isnumeric(x) && ...
    isscalar(x) && x > 0);
addParameter(parser, 'ProjectionRadius', 35, @(x) isnumeric(x) && ...
    isscalar(x) && isfinite(x) && x > 0);
addParameter(parser, 'ProjectionMode', '3d', @(x) ischar(x) || isstring(x));
addParameter(parser, 'BrainTemplate', 'brainnet', @(x) ischar(x) || isstring(x));
addParameter(parser, 'SurfaceFile', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ShowTitle', false, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'OutputPng', '', @(x) ischar(x) || isstring(x));
parse(parser, sourceFile, varargin{:});

sourceFile = char(parser.Results.sourceFile);
outputPng = char(parser.Results.OutputPng);

result = load(sourceFile, 'sourceERP', 'inverse');
if ~isfield(result, 'sourceERP') || ~isfield(result, 'inverse')
    error('%s:BadSourceFile', mfilename, ...
        'sourceFile must contain sourceERP and inverse.');
end
if ~isfield(result.inverse, 'sourceXYZ') || isempty(result.inverse.sourceXYZ)
    error('%s:MissingSourceXYZ', mfilename, ...
        'The inverse struct must contain sourceXYZ.');
end

sourceERP = squeeze(result.sourceERP);
sourceXYZ = result.inverse.sourceXYZ;
if ndims(sourceERP) > 2
    sourceERP = reshape(sourceERP, size(sourceERP, 1), []);
end
if size(sourceERP, 1) ~= size(sourceXYZ, 1)
    error('%s:SourceMismatch', mfilename, ...
        'sourceERP rows (%d) do not match sourceXYZ rows (%d).', ...
        size(sourceERP, 1), size(sourceXYZ, 1));
end

timeIndex = parser.Results.TimeIndex;
if isempty(timeIndex)
    [~, timeIndex] = max(mean(abs(sourceERP), 1));
else
    timeIndex = round(timeIndex);
end
if timeIndex > size(sourceERP, 2)
    error('%s:BadTimeIndex', mfilename, ...
        'TimeIndex %d exceeds available samples %d.', timeIndex, size(sourceERP, 2));
end

[vertices, faces, surfaceKind] = load_cortex_surface(char(parser.Results.BrainTemplate), ...
    char(parser.Results.SurfaceFile));

sourceValues = abs(sourceERP(:, timeIndex));
regionMask = source_region_mask(sourceXYZ, parser.Results.Region, parser.Results.ACCBox);
sourceValues(~regionMask) = 0;
sourceValues = keep_strongest_sources(sourceValues, parser.Results.MaxSources);
vertexValues = project_active_sources_to_vertices(vertices, sourceXYZ, sourceValues, ...
    parser.Results.ProjectionRadius, parser.Results.ProjectionMode);
vertexValues = smooth_vertex_values(faces, vertexValues, 3);

activeVertexValues = vertexValues(vertexValues > 0);
if isempty(activeVertexValues)
    error('%s:EmptyRegion', mfilename, ...
        'No nonzero values remain after applying Region=%s.', char(parser.Results.Region));
end

threshold = local_percentile(activeVertexValues, parser.Results.ThresholdPercentile);
colorMax = local_percentile(activeVertexValues, parser.Results.ColorPercentile);
if colorMax <= threshold
    colorMax = max(vertexValues);
end
if colorMax <= threshold
    colorMax = threshold + eps;
end

colors = cortex_colors(vertexValues, threshold, colorMax);

viewSpec = parser.Results.View;
fig = figure('Color', 'k', 'Name', sprintf('sLORETA cortex sample %d', timeIndex));
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
        annotation(fig, 'textbox', [0.02 0.955 0.96 0.04], ...
            'String', sprintf('sLORETA cortex map, time sample %d', timeIndex), ...
            'Color', 'w', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
            'FontSize', 14);
    end
else
    set(fig, 'Position', [100 100 1100 820]);
    ax = axes('Parent', fig, 'Color', 'k');
    set(ax, 'Position', [0 0 1 1]);
    render_cortex_axes(ax, vertices, faces, colors, viewSpec, 1.12, surfaceKind);
    if parser.Results.ShowTitle
        title(ax, sprintf('sLORETA cortex map, time sample %d', timeIndex), ...
            'Color', 'w', 'FontWeight', 'normal');
    end
end

if ~isempty(outputPng)
    save_figure_png(fig, outputPng, 250);
end

fprintf('Rendered cortical source map at time sample %d.\n', timeIndex);

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

function sourceValues = keep_strongest_sources(sourceValues, maxSources)
if isinf(maxSources)
    return;
end

idx = find(sourceValues > 0);
if numel(idx) <= maxSources
    return;
end

[~, order] = sort(sourceValues(idx), 'descend');
keep = idx(order(1:round(maxSources)));
maskedValues = zeros(size(sourceValues));
maskedValues(keep) = sourceValues(keep);
sourceValues = maskedValues;
end

function mask = source_region_mask(sourceXYZ, region, accBox)
switch lower(char(region))
    case 'all'
        mask = true(size(sourceXYZ, 1), 1);
    case 'acc'
        accBox = reshape(accBox, 1, 6);
        mask = sourceXYZ(:, 1) >= accBox(1) & sourceXYZ(:, 1) <= accBox(2) & ...
            sourceXYZ(:, 2) >= accBox(3) & sourceXYZ(:, 2) <= accBox(4) & ...
            sourceXYZ(:, 3) >= accBox(5) & sourceXYZ(:, 3) <= accBox(6);
    otherwise
        error('%s:BadRegion', mfilename, ...
            'Region must be ''all'' or ''acc''.');
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

function vertexValues = project_active_sources_to_vertices(vertices, sourceXYZ, sourceValues, radius, mode)
activeIdx = find(sourceValues > 0);
if isempty(activeIdx)
    vertexValues = zeros(size(vertices, 1), 1);
    return;
end

activeXYZ = sourceXYZ(activeIdx, :);
activeValues = sourceValues(activeIdx);
nVertices = size(vertices, 1);
vertexValues = zeros(nVertices, 1);
chunkSize = 1500;
sigma2 = 2 * radius ^ 2;

mode = lower(char(mode));
switch mode
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
end

if max(vertexValues) > 0
    vertexValues = vertexValues ./ max(vertexValues) .* max(activeValues);
end
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

function colors = cortex_colors(values, threshold, colorMax)
baseGray = repmat([0.58 0.58 0.56], numel(values), 1);
hotMap = hot(256);
scaled = (values - threshold) ./ (colorMax - threshold);
scaled = min(max(scaled, 0), 1);
colorIdx = max(1, min(256, round(1 + scaled * 255)));
hotColors = hotMap(colorIdx, :);
blend = scaled .^ 0.6;
colors = baseGray .* (1 - blend) + hotColors .* blend;
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
                'View must be ''left'', ''top'', ''right'', ''front'', ''back'', ''bottom'', ''six'', or [az el].');
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
                'View must be ''left'', ''top'', ''right'', ''front'', ''back'', ''bottom'', ''six'', or [az el].');
    end
end
end

function value = local_percentile(x, pct)
x = sort(x(:));
if isempty(x)
    value = NaN;
    return;
end
position = 1 + (numel(x) - 1) * pct / 100;
lowerIdx = floor(position);
upperIdx = ceil(position);
if lowerIdx == upperIdx
    value = x(lowerIdx);
else
    weight = position - lowerIdx;
    value = (1 - weight) * x(lowerIdx) + weight * x(upperIdx);
end
end
