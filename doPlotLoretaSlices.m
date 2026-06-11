function outputs = doPlotLoretaSlices(inputMatFile, varargin)
%DOPLOTLORETTASLICES Plot orthogonal internal sLORETA source slices.
%
%   outputs = doPlotLoretaSlices(inputMatFile)
%   outputs = doPlotLoretaSlices(inputMatFile, 'SliceXYZ', [0 0 35])
%   outputs = doPlotLoretaSlices(inputMatFile, 'SliceX', -10, 'SliceY', 20)
%
%   inputMatFile can be a subject-level output from doPlotsLoreta, a
%   subject-level contrast output from doContrastLoreta, or a group contrast
%   output from doGroupContrastLoreta.

parser = inputParser;
parser.FunctionName = mfilename;
addRequired(parser, 'inputMatFile', @(x) ischar(x) || isstring(x));
addParameter(parser, 'MapType', 'auto', @(x) ischar(x) || isstring(x));
addParameter(parser, 'TimeIndex', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1));
addParameter(parser, 'SliceXYZ', [], @(x) isempty(x) || ...
    (isnumeric(x) && numel(x) == 3 && all(isfinite(x))));
addParameter(parser, 'SliceX', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x)));
addParameter(parser, 'SliceY', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x)));
addParameter(parser, 'SliceZ', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x)));
addParameter(parser, 'GridStep', 3, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x > 0);
addParameter(parser, 'SlabMm', 8, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x > 0);
addParameter(parser, 'ThresholdPercentile', 70, @(x) isnumeric(x) && ...
    isscalar(x) && isfinite(x) && x >= 0 && x <= 100);
addParameter(parser, 'ColorPercentile', 99, @(x) isnumeric(x) && ...
    isscalar(x) && isfinite(x) && x >= 0 && x <= 100);
addParameter(parser, 'OutputPng', '', @(x) ischar(x) || isstring(x));
parse(parser, inputMatFile, varargin{:});
opts = parser.Results;

[sourceXYZ, values, mapLabel, isSigned, timeIndex] = load_loreta_map( ...
    char(inputMatFile), char(opts.MapType), opts.TimeIndex);
sliceXYZ = choose_slice_xyz(sourceXYZ, values, isSigned, opts);

fig = figure('Color', 'w', 'Name', sprintf('sLORETA slices: %s', mapLabel));
set(fig, 'Position', [100 100 1500 500]);

plot_slice_panel(fig, 1, 'sagittal', sourceXYZ, values, sliceXYZ, opts, isSigned);
plot_slice_panel(fig, 2, 'coronal', sourceXYZ, values, sliceXYZ, opts, isSigned);
plot_slice_panel(fig, 3, 'axial', sourceXYZ, values, sliceXYZ, opts, isSigned);

if ~isempty(char(opts.OutputPng))
    save_figure_png(fig, char(opts.OutputPng), 250);
end

outputs = struct();
outputs.figure = fig;
outputs.outputPng = char(opts.OutputPng);
outputs.inputMatFile = char(inputMatFile);
outputs.mapLabel = mapLabel;
outputs.isSigned = isSigned;
outputs.timeIndex = timeIndex;
outputs.sliceXYZ = sliceXYZ;

fprintf('sLORETA slice plot rendered through XYZ [%s].\n', num2str(sliceXYZ));
if ~isempty(char(opts.OutputPng))
    fprintf('sLORETA slice plot saved: %s\n', char(opts.OutputPng));
end

end

function [sourceXYZ, values, mapLabel, isSigned, timeIndex] = load_loreta_map(inputMatFile, mapType, timeIndex)
loaded = load(inputMatFile);
mapType = lower(mapType);
timeIndex = select_time_index(loaded, timeIndex);

if strcmp(mapType, 'auto')
    if isfield(loaded, 'groupSourceStats')
        mapType = 'meancontrast';
    elseif isfield(loaded, 'sourceContrast')
        mapType = 'contrast';
    else
        mapType = 'activity';
    end
end

switch mapType
    case {'activity', 'source', 'norm'}
        require_fields(loaded, {'sourceValues', 'inverse'}, inputMatFile);
        values = loaded.sourceValues(:);
        sourceXYZ = loaded.inverse.sourceXYZ;
        mapLabel = 'activity';
        isSigned = false;
    case {'contrast', 'sourcecontrast'}
        require_fields(loaded, {'sourceValues', 'inverse'}, inputMatFile);
        values = loaded.sourceValues(:);
        sourceXYZ = loaded.inverse.sourceXYZ;
        mapLabel = 'condition 1 - condition 2';
        isSigned = true;
    case {'meancontrast', 'groupmean', 'group'}
        require_fields(loaded, {'groupSourceStats', 'inverse'}, inputMatFile);
        values = loaded.groupSourceStats.mean(:);
        sourceXYZ = loaded.inverse.sourceXYZ;
        mapLabel = 'group mean contrast';
        isSigned = true;
    case {'t', 'tvalue', 'tstat'}
        require_fields(loaded, {'groupSourceStats', 'inverse'}, inputMatFile);
        values = loaded.groupSourceStats.t(:);
        sourceXYZ = loaded.inverse.sourceXYZ;
        mapLabel = 'group t value';
        isSigned = true;
    case {'p', 'pmap', 'signedlogp'}
        require_fields(loaded, {'groupSourceStats', 'inverse'}, inputMatFile);
        values = sign(loaded.groupSourceStats.mean(:)) .* ...
            -log10(max(loaded.groupSourceStats.p(:), realmin));
        sourceXYZ = loaded.inverse.sourceXYZ;
        mapLabel = 'signed -log10(p)';
        isSigned = true;
    case {'significant', 'sig'}
        require_fields(loaded, {'groupSourceStats', 'inverse'}, inputMatFile);
        values = loaded.groupSourceStats.t(:);
        values(~loaded.groupSourceStats.significant(:)) = 0;
        sourceXYZ = loaded.inverse.sourceXYZ;
        mapLabel = 'significant group t value';
        isSigned = true;
    otherwise
        error([mfilename ':BadMapType'], ...
            'Unknown MapType %s.', mapType);
end
end

function timeIndex = select_time_index(loaded, timeIndex)
if ~isempty(timeIndex)
    timeIndex = round(timeIndex);
    return;
end
if isfield(loaded, 'metadata') && isfield(loaded.metadata, 'timeIndex')
    timeIndex = loaded.metadata.timeIndex;
else
    timeIndex = [];
end
end

function require_fields(loaded, fieldNames, inputMatFile)
for idx = 1:numel(fieldNames)
    if ~isfield(loaded, fieldNames{idx})
        error([mfilename ':MissingField'], ...
            '%s does not contain required field %s.', inputMatFile, fieldNames{idx});
    end
end
end

function sliceXYZ = choose_slice_xyz(sourceXYZ, values, isSigned, opts)
if isempty(opts.SliceXYZ)
    if isSigned
        [~, peakIdx] = max(abs(values));
    else
        [~, peakIdx] = max(values);
    end
    sliceXYZ = sourceXYZ(peakIdx, :);
else
    sliceXYZ = opts.SliceXYZ(:).';
end
if ~isempty(opts.SliceX)
    sliceXYZ(1) = opts.SliceX;
end
if ~isempty(opts.SliceY)
    sliceXYZ(2) = opts.SliceY;
end
if ~isempty(opts.SliceZ)
    sliceXYZ(3) = opts.SliceZ;
end
sliceXYZ = clamp_slice_xyz(sliceXYZ, sourceXYZ);
end

function sliceXYZ = clamp_slice_xyz(sliceXYZ, sourceXYZ)
mins = min(sourceXYZ, [], 1);
maxs = max(sourceXYZ, [], 1);
sliceXYZ = max(mins, min(maxs, sliceXYZ));
end

function plot_slice_panel(fig, panelIdx, planeName, sourceXYZ, values, sliceXYZ, opts, isSigned)
ax = subplot(1, 3, panelIdx, 'Parent', fig);
[gridA, gridB, overlay, axisLabels, titleText] = slice_overlay( ...
    planeName, sourceXYZ, values, sliceXYZ, opts.GridStep, opts.SlabMm);

imagesc(ax, gridA(1, :), gridB(:, 1), overlay);
axis(ax, 'image');
axis(ax, 'xy');
box(ax, 'off');
xlabel(ax, axisLabels{1});
ylabel(ax, axisLabels{2});
title(ax, titleText, 'FontWeight', 'normal');
hold(ax, 'on');
plot_crosshair(ax, planeName, sliceXYZ);
style_axis(ax);
apply_colors(ax, values, overlay, opts, isSigned);
end

function [gridA, gridB, overlay, axisLabels, titleText] = slice_overlay( ...
    planeName, sourceXYZ, values, sliceXYZ, gridStep, slabMm)
switch planeName
    case 'sagittal'
        distance = abs(sourceXYZ(:, 1) - sliceXYZ(1));
        pointsA = sourceXYZ(:, 2);
        pointsB = sourceXYZ(:, 3);
        rangeA = min(sourceXYZ(:, 2)):gridStep:max(sourceXYZ(:, 2));
        rangeB = min(sourceXYZ(:, 3)):gridStep:max(sourceXYZ(:, 3));
        axisLabels = {'y (mm)', 'z (mm)'};
        titleText = sprintf('Sagittal x = %.1f mm', sliceXYZ(1));
    case 'coronal'
        distance = abs(sourceXYZ(:, 2) - sliceXYZ(2));
        pointsA = sourceXYZ(:, 1);
        pointsB = sourceXYZ(:, 3);
        rangeA = min(sourceXYZ(:, 1)):gridStep:max(sourceXYZ(:, 1));
        rangeB = min(sourceXYZ(:, 3)):gridStep:max(sourceXYZ(:, 3));
        axisLabels = {'x (mm)', 'z (mm)'};
        titleText = sprintf('Coronal y = %.1f mm', sliceXYZ(2));
    case 'axial'
        distance = abs(sourceXYZ(:, 3) - sliceXYZ(3));
        pointsA = sourceXYZ(:, 1);
        pointsB = sourceXYZ(:, 2);
        rangeA = min(sourceXYZ(:, 1)):gridStep:max(sourceXYZ(:, 1));
        rangeB = min(sourceXYZ(:, 2)):gridStep:max(sourceXYZ(:, 2));
        axisLabels = {'x (mm)', 'y (mm)'};
        titleText = sprintf('Axial z = %.1f mm', sliceXYZ(3));
end

[gridA, gridB] = meshgrid(rangeA, rangeB);
keep = distance <= slabMm;
if nnz(keep) < 6
    keep = distance <= slabMm * 2;
end
if nnz(keep) < 3
    [~, order] = sort(distance, 'ascend');
    keep(order(1:min(12, numel(order)))) = true;
end

weights = exp(-(distance(keep) .^ 2) ./ (2 * max(slabMm / 2, eps) ^ 2));
weightedValues = values(keep) .* weights;
overlay = griddata(pointsA(keep), pointsB(keep), weightedValues, ...
    gridA, gridB, 'linear');
nearestOverlay = griddata(pointsA(keep), pointsB(keep), weightedValues, ...
    gridA, gridB, 'nearest');
overlay(isnan(overlay)) = nearestOverlay(isnan(overlay));
overlay(isnan(overlay)) = 0;
end

function plot_crosshair(ax, planeName, sliceXYZ)
yl = ylim(ax);
xl = xlim(ax);
switch planeName
    case 'sagittal'
        plot(ax, [sliceXYZ(2) sliceXYZ(2)], yl, 'k-', 'LineWidth', 0.5);
        plot(ax, xl, [sliceXYZ(3) sliceXYZ(3)], 'k-', 'LineWidth', 0.5);
    case 'coronal'
        plot(ax, [sliceXYZ(1) sliceXYZ(1)], yl, 'k-', 'LineWidth', 0.5);
        plot(ax, xl, [sliceXYZ(3) sliceXYZ(3)], 'k-', 'LineWidth', 0.5);
    case 'axial'
        plot(ax, [sliceXYZ(1) sliceXYZ(1)], yl, 'k-', 'LineWidth', 0.5);
        plot(ax, xl, [sliceXYZ(2) sliceXYZ(2)], 'k-', 'LineWidth', 0.5);
end
end

function style_axis(ax)
set(ax, 'Color', [0.96 0.96 0.96], ...
    'XColor', [0.2 0.2 0.2], ...
    'YColor', [0.2 0.2 0.2], ...
    'FontSize', 10);
end

function apply_colors(ax, values, overlay, opts, isSigned)
if isSigned
    limit = local_percentile(abs(values), opts.ColorPercentile);
    if limit <= 0
        limit = max(abs(overlay(:)));
    end
    if limit <= 0
        limit = 1;
    end
    colormap(ax, blue_white_red(256));
    caxis(ax, [-limit limit]);
else
    threshold = local_percentile(values, opts.ThresholdPercentile);
    colorMax = local_percentile(values, opts.ColorPercentile);
    if colorMax <= threshold
        colorMax = max(values);
    end
    if colorMax <= 0
        colorMax = 1;
    end
    colormap(ax, hot(256));
    caxis(ax, [0 colorMax]);
end
cb = colorbar(ax);
set(cb, 'Box', 'off');
end

function cmap = blue_white_red(n)
if nargin < 1
    n = 256;
end
half = floor(n / 2);
blue = [linspace(0, 1, half).', linspace(0.15, 1, half).', ones(half, 1)];
red = [ones(n - half, 1), linspace(1, 0.12, n - half).', linspace(1, 0, n - half).'];
cmap = [blue; red];
end

function save_figure_png(fig, outputPng, resolution)
parentDir = fileparts(outputPng);
if ~isempty(parentDir) && ~exist(parentDir, 'dir')
    mkdir(parentDir);
end
if exist('exportgraphics', 'file') == 2
    exportgraphics(fig, outputPng, 'Resolution', resolution, ...
        'BackgroundColor', 'white');
else
    set(fig, 'InvertHardcopy', 'off');
    print(fig, outputPng, '-dpng', sprintf('-r%d', resolution));
end
end

function value = local_percentile(x, pct)
x = sort(abs(x(:)));
x = x(isfinite(x));
if isempty(x)
    value = 0;
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
