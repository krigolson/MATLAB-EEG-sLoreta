function outputs = doPlotLoretaMriSlices(inputMatFile, varargin)
%DOPLOTLORETTAMRISLICES Overlay sLORETA source values on template MRI slices.
%
%   outputs = doPlotLoretaMriSlices(inputMatFile)
%   outputs = doPlotLoretaMriSlices(inputMatFile, 'SliceXYZ', [0 0 35])
%
%   inputMatFile can be a subject-level output, a subject contrast output,
%   or a group contrast output from this toolbox.

parser = inputParser;
parser.FunctionName = mfilename;
addRequired(parser, 'inputMatFile', @(x) ischar(x) || isstring(x));
addParameter(parser, 'MapType', 'auto', @(x) ischar(x) || isstring(x));
addParameter(parser, 'TimeIndex', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1));
addParameter(parser, 'SliceXYZ', [], @(x) isempty(x) || ...
    (isnumeric(x) && size(x, 2) == 3 && all(isfinite(x(:)))));
addParameter(parser, 'SliceX', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x)));
addParameter(parser, 'SliceY', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x)));
addParameter(parser, 'SliceZ', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x)));
addParameter(parser, 'SlabMm', 12, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x > 0);
addParameter(parser, 'AutoSliceMinZ', -35, @(x) isnumeric(x) && isscalar(x));
addParameter(parser, 'OverlaySigmaMm', 10, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x > 0);
addParameter(parser, 'NumSliceRows', 2, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x >= 1 && x <= 4);
addParameter(parser, 'MinPeakDistanceMm', 35, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x >= 0);
addParameter(parser, 'ThresholdPercentile', 80, @(x) isnumeric(x) && ...
    isscalar(x) && isfinite(x) && x >= 0 && x <= 100);
addParameter(parser, 'ColorPercentile', 99, @(x) isnumeric(x) && ...
    isscalar(x) && isfinite(x) && x >= 0 && x <= 100);
addParameter(parser, 'TemplateMriFile', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'TemplateSegFile', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'BrainMaskLabel', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(parser, 'BrainMaskDilateVoxels', 2, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x >= 0);
addParameter(parser, 'OutputPng', '', @(x) ischar(x) || isstring(x));
parse(parser, inputMatFile, varargin{:});
opts = parser.Results;

[sourceXYZ, sourceValues, mapLabel, isSigned, timeIndex] = load_loreta_map( ...
    char(inputMatFile), char(opts.MapType), opts.TimeIndex);
sliceXYZ = choose_slice_xyz(sourceXYZ, sourceValues, isSigned, opts);

mriFile = char(opts.TemplateMriFile);
if isempty(mriFile)
    mriFile = fullfile(fileparts(mfilename('fullpath')), ...
        'templates', 'dipfit', 'standard_BEM', 'standard_mri.mat');
else
    mriFile = resolve_template_path(mriFile);
end
if ~exist(mriFile, 'file')
    error([mfilename ':MissingMRI'], ...
        'Could not find template MRI file: %s', mriFile);
end
mriStruct = load(mriFile, 'mri');
mri = mriStruct.mri;
anatomy = double(mri.anatomy);
brainMaskVolume = load_brain_mask_volume(opts, size(anatomy));
worldToVoxel = inv(mri.transform);

if isSigned
    threshold = local_percentile(abs(sourceValues), opts.ThresholdPercentile);
    colorMax = local_percentile(abs(sourceValues), opts.ColorPercentile);
    if colorMax <= threshold
        colorMax = max(abs(sourceValues));
    end
else
    threshold = local_percentile(sourceValues, opts.ThresholdPercentile);
    colorMax = local_percentile(sourceValues, opts.ColorPercentile);
    if colorMax <= threshold
        colorMax = max(sourceValues);
    end
end
if colorMax <= threshold
    colorMax = threshold + eps;
end

numRows = size(sliceXYZ, 1);
fig = figure('Color', 'k', 'Name', sprintf('sLORETA MRI slices: %s', mapLabel));
set(fig, 'Position', [100 100 1320 430 * numRows]);

for rowIdx = 1:numRows
    sliceVoxel = round(worldToVoxel * [sliceXYZ(rowIdx, :).'; 1]);
    sliceVoxel = sliceVoxel(1:3).';
    sliceVoxel = max([1 1 1], min(size(anatomy), sliceVoxel));

    panelBase = (rowIdx - 1) * 3;
    rowLabel = sprintf('Peak %d', rowIdx);
    plot_sagittal(fig, panelBase + 1, numRows, anatomy, mri.transform, sliceVoxel, ...
        brainMaskVolume, sourceXYZ, sourceValues, opts.SlabMm, opts.OverlaySigmaMm, ...
        threshold, colorMax, isSigned, rowLabel);
    plot_coronal(fig, panelBase + 2, numRows, anatomy, mri.transform, sliceVoxel, ...
        brainMaskVolume, sourceXYZ, sourceValues, opts.SlabMm, opts.OverlaySigmaMm, ...
        threshold, colorMax, isSigned, rowLabel);
    plot_axial(fig, panelBase + 3, numRows, anatomy, mri.transform, sliceVoxel, ...
        brainMaskVolume, sourceXYZ, sourceValues, opts.SlabMm, opts.OverlaySigmaMm, ...
        threshold, colorMax, isSigned, rowLabel);
end

if ~isempty(char(opts.OutputPng))
    save_figure_png(fig, char(opts.OutputPng), 250);
end

if isempty(timeIndex)
    fprintf('MRI slice overlay rendered through XYZ rows:\n');
else
    fprintf('MRI slice overlay rendered at sample %d through XYZ rows:\n', timeIndex);
end
disp(sliceXYZ);

outputs = struct();
outputs.figure = fig;
outputs.outputPng = char(opts.OutputPng);
outputs.inputMatFile = char(inputMatFile);
outputs.mapLabel = mapLabel;
outputs.isSigned = isSigned;
outputs.timeIndex = timeIndex;
outputs.sliceXYZ = sliceXYZ;

if ~isempty(char(opts.OutputPng))
    fprintf('MRI slice overlay saved: %s\n', char(opts.OutputPng));
end

end

function plot_sagittal(fig, panelIdx, numRows, anatomy, transform, sliceVoxel, ...
    brainMaskVolume, sourceXYZ, sourceValues, slabMm, overlaySigmaMm, threshold, colorMax, isSigned, rowLabel)
anatSlice = squeeze(anatomy(sliceVoxel(1), :, :)).';
brainMask = squeeze(brainMaskVolume(sliceVoxel(1), :, :)).';
[worldY, worldZ] = plane_world_grid(size(anatSlice), transform, sliceVoxel, 'sagittal');
sliceWorld = voxel_to_world(transform, sliceVoxel, []);
overlay = plane_overlay(sourceXYZ, sourceValues, worldY, worldZ, ...
    sliceVoxel, transform, 'sagittal', slabMm, overlaySigmaMm);
[anatSlice, brainMask, overlay, worldY, worldZ] = flip_slice_rows( ...
    anatSlice, brainMask, overlay, worldY, worldZ);
plot_panel(fig, panelIdx, anatSlice, brainMask, overlay, threshold, colorMax, ...
    sprintf('%s  sagittal  x = %.1f', rowLabel, voxel_to_world(transform, sliceVoxel, 1)), ...
    isSigned, worldY, worldZ, sliceWorld(2), sliceWorld(3), numRows);
end

function plot_coronal(fig, panelIdx, numRows, anatomy, transform, sliceVoxel, ...
    brainMaskVolume, sourceXYZ, sourceValues, slabMm, overlaySigmaMm, threshold, colorMax, isSigned, rowLabel)
anatSlice = squeeze(anatomy(:, sliceVoxel(2), :)).';
brainMask = squeeze(brainMaskVolume(:, sliceVoxel(2), :)).';
[worldX, worldZ] = plane_world_grid(size(anatSlice), transform, sliceVoxel, 'coronal');
sliceWorld = voxel_to_world(transform, sliceVoxel, []);
overlay = plane_overlay(sourceXYZ, sourceValues, worldX, worldZ, ...
    sliceVoxel, transform, 'coronal', slabMm, overlaySigmaMm);
[anatSlice, brainMask, overlay, worldX, worldZ] = flip_slice_rows( ...
    anatSlice, brainMask, overlay, worldX, worldZ);
plot_panel(fig, panelIdx, anatSlice, brainMask, overlay, threshold, colorMax, ...
    sprintf('%s  coronal  y = %.1f', rowLabel, voxel_to_world(transform, sliceVoxel, 2)), ...
    isSigned, worldX, worldZ, sliceWorld(1), sliceWorld(3), numRows);
end

function plot_axial(fig, panelIdx, numRows, anatomy, transform, sliceVoxel, ...
    brainMaskVolume, sourceXYZ, sourceValues, slabMm, overlaySigmaMm, threshold, colorMax, isSigned, rowLabel)
anatSlice = squeeze(anatomy(:, :, sliceVoxel(3))).';
brainMask = squeeze(brainMaskVolume(:, :, sliceVoxel(3))).';
[worldX, worldY] = plane_world_grid(size(anatSlice), transform, sliceVoxel, 'axial');
sliceWorld = voxel_to_world(transform, sliceVoxel, []);
overlay = plane_overlay(sourceXYZ, sourceValues, worldX, worldY, ...
    sliceVoxel, transform, 'axial', slabMm, overlaySigmaMm);
[anatSlice, brainMask, overlay, worldX, worldY] = flip_slice_rows( ...
    anatSlice, brainMask, overlay, worldX, worldY);
plot_panel(fig, panelIdx, anatSlice, brainMask, overlay, threshold, colorMax, ...
    sprintf('%s  axial  z = %.1f', rowLabel, voxel_to_world(transform, sliceVoxel, 3)), ...
    isSigned, worldX, worldY, sliceWorld(1), sliceWorld(2), numRows);
end

function [anatSlice, brainMask, overlay, worldA, worldB] = flip_slice_rows( ...
    anatSlice, brainMask, overlay, worldA, worldB)
anatSlice = flipud(anatSlice);
brainMask = flipud(brainMask);
overlay = flipud(overlay);
worldA = flipud(worldA);
worldB = flipud(worldB);
end

function plot_panel(fig, panelIdx, anatSlice, brainMask, overlay, threshold, colorMax, titleText, isSigned, worldA, worldB, crossA, crossB, numRows)
ax = subplot(numRows, 3, panelIdx, 'Parent', fig);
set(ax, 'Position', compact_panel_position(panelIdx, numRows));
baseImg = anatomy_rgb(anatSlice, brainMask);
overlay(~brainMask) = 0;
visibleValues = overlay(brainMask & isfinite(overlay) & overlay ~= 0);
if ~isempty(visibleValues)
    if isSigned
        localThreshold = local_percentile(abs(visibleValues), 82);
        localColorMax = local_percentile(abs(visibleValues), 99.5);
    else
        localThreshold = local_percentile(visibleValues, 82);
        localColorMax = local_percentile(visibleValues, 99.5);
    end
    if isfinite(localThreshold) && isfinite(localColorMax) && localColorMax > localThreshold
        threshold = max(0, localThreshold);
        colorMax = localColorMax;
    end
end
if isSigned
    cmap = blue_white_red(256);
    scaled = min(max((overlay + colorMax) ./ (2 * colorMax), 0), 1);
    alphaData = min(max((abs(overlay) - threshold) ./ (colorMax - threshold), 0), 1);
else
    cmap = hot(256);
    scaled = min(max((overlay - threshold) ./ (colorMax - threshold), 0), 1);
    alphaData = scaled;
end

overlayImg = ind2rgb(uint8(1 + 255 * scaled), cmap);
alphaData = 0.68 * (alphaData .^ 0.8);
alphaData(~isfinite(alphaData)) = 0;
alphaData = min(max(alphaData, 0), 0.9);
alpha3 = repmat(alphaData, [1 1 3]);
compositeImg = baseImg .* (1 - alpha3) + overlayImg .* alpha3;

image(ax, compositeImg);
axis(ax, 'image');
axis(ax, 'off');
hold(ax, 'on');
draw_crosshair(ax, worldA, worldB, crossA, crossB, size(compositeImg));

title(ax, titleText, 'Color', 'w', 'FontWeight', 'normal', 'FontSize', 10);
end

function pos = compact_panel_position(panelIdx, numRows)
leftMargin = 0.035;
rightMargin = 0.025;
topMargin = 0.055;
bottomMargin = 0.045;
hGap = 0.025;
vGap = 0.075;

colIdx = mod(panelIdx - 1, 3) + 1;
rowIdx = floor((panelIdx - 1) / 3) + 1;
panelWidth = (1 - leftMargin - rightMargin - 2 * hGap) / 3;
panelHeight = (1 - topMargin - bottomMargin - (numRows - 1) * vGap) / numRows;
left = leftMargin + (colIdx - 1) * (panelWidth + hGap);
bottom = 1 - topMargin - rowIdx * panelHeight - (rowIdx - 1) * vGap;
pos = [left bottom panelWidth panelHeight];
end

function draw_crosshair(ax, worldA, worldB, crossA, crossB, imageSize)
cols = size(worldA, 2);
rows = size(worldA, 1);
aByCol = median(worldA, 1);
bByRow = median(worldB, 2);
[~, colIdx] = min(abs(aByCol - crossA));
[~, rowIdx] = min(abs(bByRow - crossB));
colIdx = max(1, min(cols, colIdx));
rowIdx = max(1, min(rows, rowIdx));

line(ax, [1 imageSize(2)], [rowIdx rowIdx], ...
    'Color', [0.25 0.85 0.30], 'LineWidth', 0.8);
line(ax, [colIdx colIdx], [1 imageSize(1)], ...
    'Color', [0.15 0.35 1.00], 'LineWidth', 0.8);
plot(ax, colIdx, rowIdx, '+', 'Color', [1 0.25 0.25], ...
    'MarkerSize', 7, 'LineWidth', 0.9);
end

function rgb = anatomy_rgb(anatSlice, brainMask)
anatSlice = double(anatSlice);
nonzero = anatSlice(anatSlice > 0 & brainMask);
if isempty(nonzero)
    lo = min(anatSlice(:));
    hi = max(anatSlice(:));
else
    lo = local_percentile(nonzero, 1);
    hi = local_percentile(nonzero, 99.5);
end
if hi <= lo
    hi = lo + eps;
end
scaled = min(max((anatSlice - lo) ./ (hi - lo), 0), 1);
scaled = 0.08 + 0.92 * (scaled .^ 0.65);
scaled(anatSlice <= 0 | ~brainMask) = 0;
rgb = repmat(scaled, [1 1 3]);
end

function brainMaskVolume = load_brain_mask_volume(opts, anatomySize)
segFile = char(opts.TemplateSegFile);
if isempty(segFile)
    segFile = fullfile(fileparts(mfilename('fullpath')), ...
        'templates', 'dipfit', 'standard_BEM', 'standard_seg.mat');
else
    segFile = resolve_template_path(segFile);
end

if exist(segFile, 'file')
    segStruct = load(segFile, 'mri');
    if isfield(segStruct, 'mri') && isfield(segStruct.mri, 'seg')
        brainMaskVolume = segStruct.mri.seg == opts.BrainMaskLabel;
        brainMaskVolume = dilate_mask(brainMaskVolume, round(opts.BrainMaskDilateVoxels));
        if isequal(size(brainMaskVolume), anatomySize)
            return;
        end
        warning([mfilename ':BadSegmentationSize'], ...
            'TemplateSegFile size does not match TemplateMriFile. Falling back to anatomy mask.');
    end
else
    warning([mfilename ':MissingSegmentation'], ...
        'Could not find template segmentation file: %s. Falling back to anatomy mask.', segFile);
end

brainMaskVolume = true(anatomySize);
end

function mask = dilate_mask(mask, nVoxels)
mask = logical(mask);
kernel = ones(3, 3, 3);
for idx = 1:nVoxels
    mask = convn(double(mask), kernel, 'same') > 0;
end
end

function filePath = resolve_template_path(filePath)
filePath = char(filePath);
if is_absolute_path(filePath) || exist(filePath, 'file')
    return;
end
filePath = fullfile(fileparts(mfilename('fullpath')), filePath);
end

function tf = is_absolute_path(filePath)
filePath = char(filePath);
tf = startsWith(filePath, filesep) || ...
    (~isempty(regexp(filePath, '^[A-Za-z]:[\\/]', 'once')));
end

function sliceXYZ = choose_slice_xyz(sourceXYZ, values, isSigned, opts)
if isempty(opts.SliceXYZ)
    candidate = true(size(values));
    if isempty(opts.SliceZ) && isfinite(opts.AutoSliceMinZ)
        candidate = candidate & sourceXYZ(:, 3) >= opts.AutoSliceMinZ;
    end
    if ~any(candidate)
        candidate = true(size(values));
    end
    if isSigned
        score = abs(values);
    else
        score = values;
    end
    numRows = round(opts.NumSliceRows);
    sliceXYZ = zeros(numRows, 3);
    usedIdx = [];
    for rowIdx = 1:numRows
        rowCandidate = candidate;
        for used = usedIdx(:).'
            distanceFromUsed = sqrt(sum((sourceXYZ - sourceXYZ(used, :)) .^ 2, 2));
            rowCandidate = rowCandidate & distanceFromUsed >= opts.MinPeakDistanceMm;
        end
        if ~any(rowCandidate)
            rowCandidate = candidate;
        end
        candidateIdx = find(rowCandidate);
        [~, localPeakIdx] = max(score(candidateIdx));
        peakIdx = candidateIdx(localPeakIdx);
        usedIdx(end + 1) = peakIdx; %#ok<AGROW>
        sliceXYZ(rowIdx, :) = sourceXYZ(peakIdx, :);
    end
else
    sliceXYZ = reshape(opts.SliceXYZ, [], 3);
end
if ~isempty(opts.SliceX)
    sliceXYZ(:, 1) = opts.SliceX;
end
if ~isempty(opts.SliceY)
    sliceXYZ(:, 2) = opts.SliceY;
end
if ~isempty(opts.SliceZ)
    sliceXYZ(:, 3) = opts.SliceZ;
end
mins = min(sourceXYZ, [], 1);
maxs = max(sourceXYZ, [], 1);
sliceXYZ = bsxfun(@max, mins, bsxfun(@min, maxs, sliceXYZ));
end

function mask = plane_brain_mask(brainBoundary, worldA, worldB, planeName, sliceCoord, slabMm)
if isempty(brainBoundary)
    mask = true(size(worldA));
    return;
end
switch planeName
    case 'sagittal'
        distance = abs(brainBoundary(:, 1) - sliceCoord);
        pointsA = brainBoundary(:, 2);
        pointsB = brainBoundary(:, 3);
    case 'coronal'
        distance = abs(brainBoundary(:, 2) - sliceCoord);
        pointsA = brainBoundary(:, 1);
        pointsB = brainBoundary(:, 3);
    case 'axial'
        distance = abs(brainBoundary(:, 3) - sliceCoord);
        pointsA = brainBoundary(:, 1);
        pointsB = brainBoundary(:, 2);
end

try
    keep = distance <= slabMm * 1.5;
    if nnz(keep) < 20
        keep = distance <= slabMm * 3;
    end
    if nnz(keep) < 20
        [~, order] = sort(distance, 'ascend');
        keep(order(1:min(80, numel(order)))) = true;
    end
    hullIdx = convhull(pointsA(keep), pointsB(keep));
    keptA = pointsA(keep);
    keptB = pointsB(keep);
    mask = inpolygon(worldA, worldB, keptA(hullIdx), keptB(hullIdx));
catch
    mask = true(size(worldA));
end
end

function overlay = plane_overlay(sourceXYZ, sourceValues, worldA, worldB, ...
    sliceVoxel, transform, planeName, slabMm, overlaySigmaMm)
sliceWorld = voxel_to_world(transform, sliceVoxel, []);

switch planeName
    case 'sagittal'
        distance = abs(sourceXYZ(:, 1) - sliceWorld(1));
        pointsA = sourceXYZ(:, 2);
        pointsB = sourceXYZ(:, 3);
    case 'coronal'
        distance = abs(sourceXYZ(:, 2) - sliceWorld(2));
        pointsA = sourceXYZ(:, 1);
        pointsB = sourceXYZ(:, 3);
    case 'axial'
        distance = abs(sourceXYZ(:, 3) - sliceWorld(3));
        pointsA = sourceXYZ(:, 1);
        pointsB = sourceXYZ(:, 2);
end

keep = distance <= max(slabMm * 2, overlaySigmaMm * 2);
if nnz(keep) < 12
    keep = distance <= max(slabMm * 4, overlaySigmaMm * 4);
end
if ~any(keep)
    overlay = zeros(size(worldA));
    return;
end

gridA = worldA(:);
gridB = worldB(:);
overlayFlat = zeros(numel(gridA), 1);
weightFlat = zeros(numel(gridA), 1);
keepIdx = find(keep);
slabSigma = max(slabMm / 2, eps);
chunkSize = 128;

for startIdx = 1:chunkSize:numel(keepIdx)
    ids = keepIdx(startIdx:min(startIdx + chunkSize - 1, numel(keepIdx)));
    da = bsxfun(@minus, gridA, pointsA(ids).');
    db = bsxfun(@minus, gridB, pointsB(ids).');
    inPlaneWeight = exp(-(da .^ 2 + db .^ 2) ./ (2 * overlaySigmaMm ^ 2));
    throughPlaneWeight = exp(-(distance(ids).' .^ 2) ./ (2 * slabSigma ^ 2));
    weights = bsxfun(@times, inPlaneWeight, throughPlaneWeight);
    overlayFlat = overlayFlat + weights * sourceValues(ids);
    weightFlat = weightFlat + sum(weights, 2);
end

valid = weightFlat > max(weightFlat) * 0.01;
overlayFlat(valid) = overlayFlat(valid) ./ weightFlat(valid);
overlayFlat(~valid) = 0;
overlay = reshape(overlayFlat, size(worldA));
overlay(~isfinite(overlay)) = 0;
end

function [worldA, worldB] = plane_world_grid(sliceSize, transform, sliceVoxel, planeName)
[cols, rows] = meshgrid(1:sliceSize(2), 1:sliceSize(1));

switch planeName
    case 'sagittal'
        vox = [repmat(sliceVoxel(1), numel(rows), 1), cols(:), rows(:)];
        world = transform * [vox, ones(size(vox, 1), 1)].';
        worldA = reshape(world(2, :), sliceSize);
        worldB = reshape(world(3, :), sliceSize);
    case 'coronal'
        vox = [cols(:), repmat(sliceVoxel(2), numel(rows), 1), rows(:)];
        world = transform * [vox, ones(size(vox, 1), 1)].';
        worldA = reshape(world(1, :), sliceSize);
        worldB = reshape(world(3, :), sliceSize);
    case 'axial'
        vox = [cols(:), rows(:), repmat(sliceVoxel(3), numel(rows), 1)];
        world = transform * [vox, ones(size(vox, 1), 1)].';
        worldA = reshape(world(1, :), sliceSize);
        worldB = reshape(world(2, :), sliceSize);
end
end

function world = voxel_to_world(transform, voxel, component)
world4 = transform * [voxel(:); 1];
world = world4(1:3).';
if ~isempty(component)
    world = world(component);
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
        values = abs(loaded.sourceValues(:));
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
        'BackgroundColor', 'black');
else
    set(fig, 'InvertHardcopy', 'off');
    print(fig, outputPng, '-dpng', sprintf('-r%d', resolution));
end
end

function value = local_percentile(x, pct)
x = sort(x(:));
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
