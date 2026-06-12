function outputs = doContrastLoreta(dataInput, varargin)
%DOCONTRASTLORETA Run subject-level sLORETA contrast: condition 1 - condition 2.
%
%   outputs = doContrastLoreta(dataInput)
%   outputs = doContrastLoreta(dataInput, 'Name', value, ...)
%
%   dataInput can be a numeric matrix or a .mat filename. The expected data
%   shape is:
%
%       channels x time x conditions x subjects
%
%   The contrast is always:
%
%       condition1Idx - condition2Idx
%
%   Arrange your EEG matrix so condition 1 is the positive condition and
%   condition 2 is the negative condition.

if nargin < 1 || isempty(dataInput)
    dataInput = 'exampleSubjectData.mat';
end

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'DataVariable', 'exampleData', @(x) ischar(x) || isstring(x));
addParameter(parser, 'LeadfieldFile', 'leadfield.mat', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ChanlocsFile', 'matlocs.mat', @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputDir', pwd, @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputPrefix', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'SubjectIdx', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(parser, 'Condition1Idx', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(parser, 'Condition2Idx', 2, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(parser, 'TimeIndex', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1));
addParameter(parser, 'Lambda', 0.05, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(parser, 'Reference', 'average', @(x) ischar(x) || isstring(x));
addParameter(parser, 'AtlasName', 'DKT', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ProjectionRadius', 35, @(x) isnumeric(x) && ...
    isscalar(x) && isfinite(x) && x > 0);
addParameter(parser, 'ProjectionMode', '3d', @(x) ischar(x) || isstring(x));
addParameter(parser, 'BrainTemplate', 'brainnet', @(x) ischar(x) || isstring(x));
addParameter(parser, 'SurfaceFile', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ColorPercentile', 99, @(x) isnumeric(x) && ...
    isscalar(x) && x >= 0 && x <= 100);
addParameter(parser, 'MakeSlicePlot', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'MakeMriSlicePlot', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'CloseFigures', false, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'SliceXYZ', [], @(x) isempty(x) || ...
    (isnumeric(x) && size(x, 2) == 3 && all(isfinite(x(:)))));
addParameter(parser, 'SliceX', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x)));
addParameter(parser, 'SliceY', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x)));
addParameter(parser, 'SliceZ', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x)));
addParameter(parser, 'AutoSliceMinZ', -35, @(x) isnumeric(x) && isscalar(x));
addParameter(parser, 'OverlaySigmaMm', 10, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x > 0);
addParameter(parser, 'NumSliceRows', 2, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x >= 1 && x <= 4);
addParameter(parser, 'MinPeakDistanceMm', 35, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x >= 0);
addParameter(parser, 'TemplateMriFile', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'TemplateSegFile', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'BrainMaskLabel', 3, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(parser, 'BrainMaskDilateVoxels', 2, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x >= 0);
parse(parser, varargin{:});

opts = parser.Results;
rootDir = fileparts(mfilename('fullpath'));
outputDir = prepare_output_dir(char(opts.OutputDir), rootDir);

[condition1Data, condition2Data, dataInfo] = load_or_select_conditions(dataInput, opts);

leadfieldPath = resolve_path(char(opts.LeadfieldFile), rootDir);
model = load(leadfieldPath, 'leadfield', 'sourceXYZ', 'channelLabels');
if ~isfield(model, 'leadfield') || ~isfield(model, 'sourceXYZ')
    error([mfilename ':BadLeadfieldFile'], ...
        'LeadfieldFile must contain leadfield and sourceXYZ.');
end
chanlocsPath = resolve_path(char(opts.ChanlocsFile), rootDir);
channelLabels = {};
if exist(chanlocsPath, 'file')
    chan = load(chanlocsPath, 'chanlocs');
    if isfield(chan, 'chanlocs')
        channelLabels = sloreta_chanlocs_labels(chan.chanlocs);
    end
end
validate_channel_match(condition1Data, model, channelLabels, leadfieldPath);

inverse = sloreta_make_inverse(model.leadfield, ...
    'SourceXYZ', model.sourceXYZ, ...
    'Lambda', opts.Lambda, ...
    'Reference', char(opts.Reference), ...
    'Orientation', 'free', ...
    'Standardize', true);

sourceCondition1 = sloreta_apply(inverse, condition1Data, 'Output', 'norm');
sourceCondition2 = sloreta_apply(inverse, condition2Data, 'Output', 'norm');
sourceContrast = sourceCondition1 - sourceCondition2;
sourceContrast = reshape(sourceContrast, inverse.nSources, size(condition1Data, 2));

timeIndex = opts.TimeIndex;
if isempty(timeIndex)
    [~, timeIndex] = max(mean(abs(sourceContrast), 1));
else
    timeIndex = round(timeIndex);
end

sourceValues = sourceContrast(:, timeIndex);
regionContrast = rank_contrast_regions(sourceValues, inverse.sourceXYZ, opts);

if isempty(char(opts.OutputPrefix))
    outputPrefix = sprintf('sLoreta_sub%02d_cond%02dminus%02d_t%03d', ...
        opts.SubjectIdx, opts.Condition1Idx, opts.Condition2Idx, timeIndex);
else
    outputPrefix = char(opts.OutputPrefix);
end

matFile = fullfile(outputDir, [outputPrefix '_contrast.mat']);
imageFile = fullfile(outputDir, [outputPrefix '_contrast_sixview.png']);
sliceFile = fullfile(outputDir, [outputPrefix '_contrast_slices.png']);
mriSliceFile = fullfile(outputDir, [outputPrefix '_contrast_mri_slices.png']);
csvFile = fullfile(outputDir, [outputPrefix '_contrast_regions.csv']);

metadata = struct();
metadata.createdBy = mfilename;
metadata.dataInfo = dataInfo;
metadata.subjectIdx = opts.SubjectIdx;
metadata.condition1Idx = opts.Condition1Idx;
metadata.condition2Idx = opts.Condition2Idx;
metadata.contrast = 'condition1_minus_condition2';
metadata.timeIndex = timeIndex;
metadata.lambda = opts.Lambda;
metadata.reference = char(opts.Reference);
metadata.atlasName = char(opts.AtlasName);
metadata.projectionRadius = opts.ProjectionRadius;
metadata.projectionMode = char(opts.ProjectionMode);
metadata.brainTemplate = char(opts.BrainTemplate);
metadata.surfaceFile = char(opts.SurfaceFile);
metadata.autoSliceMinZ = opts.AutoSliceMinZ;
metadata.overlaySigmaMm = opts.OverlaySigmaMm;
metadata.numSliceRows = opts.NumSliceRows;
metadata.minPeakDistanceMm = opts.MinPeakDistanceMm;
metadata.brainMaskLabel = opts.BrainMaskLabel;
metadata.brainMaskDilateVoxels = opts.BrainMaskDilateVoxels;
metadata.outputPrefix = outputPrefix;

write_contrast_region_csv(csvFile, regionContrast);
save(matFile, 'sourceContrast', 'sourceCondition1', 'sourceCondition2', ...
    'sourceValues', 'inverse', 'channelLabels', 'metadata', 'regionContrast', '-v7');

colorLimit = local_percentile(abs(sourceValues), opts.ColorPercentile);
plot_sloreta_signed_cortex(inverse.sourceXYZ, sourceValues, ...
    'View', 'six', ...
    'ColorLimit', colorLimit, ...
    'ProjectionRadius', opts.ProjectionRadius, ...
    'ProjectionMode', char(opts.ProjectionMode), ...
    'BrainTemplate', opts.BrainTemplate, ...
    'SurfaceFile', opts.SurfaceFile, ...
    'OutputPng', imageFile);

if logical(opts.MakeSlicePlot)
    sliceFig = doPlotLoretaSlices(matFile, ...
        'MapType', 'contrast', ...
        'SliceXYZ', opts.SliceXYZ, ...
        'SliceX', opts.SliceX, ...
        'SliceY', opts.SliceY, ...
        'SliceZ', opts.SliceZ, ...
        'ColorPercentile', opts.ColorPercentile, ...
        'OutputPng', sliceFile);
    if logical(opts.CloseFigures)
        close(sliceFig.figure);
    end
end

if logical(opts.MakeMriSlicePlot)
    mriSliceFig = doPlotLoretaMriSlices(matFile, ...
        'MapType', 'contrast', ...
        'SliceXYZ', opts.SliceXYZ, ...
        'SliceX', opts.SliceX, ...
        'SliceY', opts.SliceY, ...
        'SliceZ', opts.SliceZ, ...
        'AutoSliceMinZ', opts.AutoSliceMinZ, ...
        'OverlaySigmaMm', opts.OverlaySigmaMm, ...
        'NumSliceRows', opts.NumSliceRows, ...
        'MinPeakDistanceMm', opts.MinPeakDistanceMm, ...
        'TemplateMriFile', opts.TemplateMriFile, ...
        'TemplateSegFile', opts.TemplateSegFile, ...
        'BrainMaskLabel', opts.BrainMaskLabel, ...
        'BrainMaskDilateVoxels', opts.BrainMaskDilateVoxels, ...
        'ColorPercentile', opts.ColorPercentile, ...
        'OutputPng', mriSliceFile);
    if logical(opts.CloseFigures)
        close(mriSliceFig.figure);
    end
end

outputs = struct();
outputs.imageFile = imageFile;
outputs.sliceFile = sliceFile;
outputs.mriSliceFile = mriSliceFile;
outputs.csvFile = csvFile;
outputs.matFile = matFile;
outputs.timeIndex = timeIndex;
outputs.sourceContrast = sourceContrast;
outputs.sourceXYZ = inverse.sourceXYZ;
outputs.regionContrast = regionContrast;
outputs.metadata = metadata;

fprintf('sLORETA contrast plot saved: %s\n', imageFile);
if logical(opts.MakeSlicePlot)
    fprintf('sLORETA contrast slice plot saved: %s\n', sliceFile);
end
if logical(opts.MakeMriSlicePlot)
    fprintf('sLORETA contrast MRI slice plot saved: %s\n', mriSliceFile);
end
fprintf('Contrast region CSV saved: %s\n', csvFile);
fprintf('Contrast MAT saved: %s\n', matFile);

end

function [condition1Data, condition2Data, info] = load_or_select_conditions(dataInput, opts)
info = struct();
if isnumeric(dataInput)
    raw = dataInput;
    info.source = 'numeric input';
else
    dataFile = char(dataInput);
    varName = char(opts.DataVariable);
    loaded = load(dataFile, varName);
    raw = loaded.(varName);
    info.source = dataFile;
    info.variable = varName;
end

if ndims(raw) < 3
    error([mfilename ':NeedConditions'], ...
        'Contrast analysis requires a conditions dimension.');
end

if ndims(raw) == 3
    condition1Data = raw(:, :, opts.Condition1Idx);
    condition2Data = raw(:, :, opts.Condition2Idx);
else
    condition1Data = raw(:, :, opts.Condition1Idx, opts.SubjectIdx);
    condition2Data = raw(:, :, opts.Condition2Idx, opts.SubjectIdx);
end

condition1Data = reshape(condition1Data, size(raw, 1), size(raw, 2), 1, 1);
condition2Data = reshape(condition2Data, size(raw, 1), size(raw, 2), 1, 1);
info.rawSize = size(raw);
info.condition1Size = size(condition1Data);
info.condition2Size = size(condition2Data);
end

function pathOut = resolve_path(pathIn, rootDir)
pathIn = char(pathIn);
if is_absolute_path(pathIn) || exist(pathIn, 'file') || exist(pathIn, 'dir')
    pathOut = pathIn;
else
    pathOut = fullfile(rootDir, pathIn);
end
end

function outputDir = prepare_output_dir(outputDir, rootDir)
outputDir = char(outputDir);
if ~is_absolute_path(outputDir)
    outputDir = fullfile(rootDir, outputDir);
end
if ~exist(outputDir, 'dir')
    [status, message] = mkdir(outputDir);
    if status ~= 1
        error([mfilename ':CannotCreateOutputDir'], ...
            'Could not create output directory %s. %s', outputDir, message);
    end
end
end

function tf = is_absolute_path(pathIn)
pathIn = char(pathIn);
tf = startsWith(pathIn, filesep) || ~isempty(regexp(pathIn, '^[A-Za-z]:[\\/]', 'once'));
end

function validate_channel_match(data, model, chanlocsLabels, leadfieldPath)
nDataChannels = size(data, 1);
nLeadfieldChannels = size(model.leadfield, 1);
if nDataChannels ~= nLeadfieldChannels
    error([mfilename ':ChannelCountMismatch'], ...
        ['Your data has %d channels, but %s expects %d channels. ', ...
        'Build or choose a leadfield file for this exact channel set and order. ', ...
        'For example: doMakeLeadfield(''ChanlocsFile'', ''matlocs_%dchan.mat'', ', ...
        '''OutputFile'', ''leadfield_%dchan.mat'').'], ...
        nDataChannels, leadfieldPath, nLeadfieldChannels, ...
        nDataChannels, nDataChannels);
end

if isfield(model, 'channelLabels') && ~isempty(model.channelLabels) && ~isempty(chanlocsLabels)
    leadfieldLabels = cellstr(model.channelLabels(:));
    chanlocsLabels = cellstr(chanlocsLabels(:));
    if numel(leadfieldLabels) == numel(chanlocsLabels) && ...
            ~all(strcmpi(leadfieldLabels, chanlocsLabels))
        firstMismatch = find(~strcmpi(leadfieldLabels, chanlocsLabels), 1);
        error([mfilename ':ChannelOrderMismatch'], ...
            ['The channel labels in ChanlocsFile do not match the channel order ', ...
            'saved in the leadfield. First mismatch at row %d: data/chanlocs=%s, ', ...
            'leadfield=%s. Use the matching ChanlocsFile and LeadfieldFile.'], ...
            firstMismatch, chanlocsLabels{firstMismatch}, ...
            leadfieldLabels{firstMismatch});
    end
end
end

function rows = rank_contrast_regions(sourceValues, sourceXYZ, opts)
cortexFile = fullfile(fileparts(mfilename('fullpath')), ...
    'templates', 'cortex', 'brainstorm_icbm152_cortex_pial_low.mat');
cortex = load(cortexFile, 'Vertices', 'Faces', 'Atlas');
vertices = cortex.Vertices;
if max(abs(vertices(:))) < 1
    vertices = vertices * 1000;
end

vertexValues = project_signed_sources_to_vertices(vertices, sourceXYZ, sourceValues, ...
    opts.ProjectionRadius, char(opts.ProjectionMode));
vertexValues = smooth_vertex_values(cortex.Faces, vertexValues, 3);

atlasIdx = find(strcmpi({cortex.Atlas.Name}, char(opts.AtlasName)), 1);
if isempty(atlasIdx)
    error([mfilename ':BadAtlas'], 'Atlas %s not found.', char(opts.AtlasName));
end

scouts = cortex.Atlas(atlasIdx).Scouts;
rows = repmat(struct('rank', [], 'atlas', '', 'label', '', 'region', '', ...
    'meanContrast', [], 'maxContrast', [], 'minContrast', [], 'absMeanContrast', [], ...
    'sumContrast', [], 'nVertices', [], 'direction', ''), numel(scouts), 1);

for idx = 1:numel(scouts)
    scoutVertices = scouts(idx).Vertices(:);
    scoutVertices = scoutVertices(scoutVertices >= 1 & scoutVertices <= numel(vertexValues));
    vals = vertexValues(scoutVertices);
    meanContrast = mean(vals);
    rows(idx).atlas = cortex.Atlas(atlasIdx).Name;
    rows(idx).label = scouts(idx).Label;
    rows(idx).region = scouts(idx).Region;
    rows(idx).meanContrast = meanContrast;
    rows(idx).maxContrast = max(vals);
    rows(idx).minContrast = min(vals);
    rows(idx).absMeanContrast = abs(meanContrast);
    rows(idx).sumContrast = sum(vals);
    rows(idx).nVertices = numel(vals);
    rows(idx).direction = contrast_direction(meanContrast);
end

[~, order] = sort([rows.absMeanContrast], 'descend');
rows = rows(order);
for idx = 1:numel(rows)
    rows(idx).rank = idx;
end
end

function write_contrast_region_csv(csvFile, rows)
ensure_parent_dir(csvFile);
fid = fopen(csvFile, 'w');
if fid < 0
    error([mfilename ':CannotWriteCSV'], 'Could not write %s.', csvFile);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, ['rank,atlas,label,region,meanContrast,maxContrast,minContrast,', ...
    'absMeanContrast,sumContrast,nVertices,direction\n']);
for idx = 1:numel(rows)
    fprintf(fid, '%d,%s,%s,%s,%.12g,%.12g,%.12g,%.12g,%.12g,%d,%s\n', ...
        rows(idx).rank, csv_escape(rows(idx).atlas), csv_escape(rows(idx).label), ...
        csv_escape(rows(idx).region), rows(idx).meanContrast, rows(idx).maxContrast, ...
        rows(idx).minContrast, rows(idx).absMeanContrast, rows(idx).sumContrast, ...
        rows(idx).nVertices, csv_escape(rows(idx).direction));
end
end

function ensure_parent_dir(filePath)
parentDir = fileparts(filePath);
if isempty(parentDir) || exist(parentDir, 'dir')
    return;
end
[status, message] = mkdir(parentDir);
if status ~= 1
    error([mfilename ':CannotCreateOutputDir'], ...
        'Could not create output directory %s. %s', parentDir, message);
end
end

function direction = contrast_direction(value)
if value > 0
    direction = 'condition1 > condition2';
elseif value < 0
    direction = 'condition2 > condition1';
else
    direction = 'no difference';
end
end

function text = csv_escape(value)
text = char(value);
text = strrep(text, '"', '""');
text = ['"', text, '"'];
end

function vertexValues = project_signed_sources_to_vertices(vertices, sourceXYZ, sourceValues, radius, mode)
activeIdx = find(sourceValues ~= 0);
activeXYZ = sourceXYZ(activeIdx, :);
activeValues = sourceValues(activeIdx);
vertexValues = zeros(size(vertices, 1), 1);
weightSums = zeros(size(vertices, 1), 1);
chunkSize = 1500;
sigma2 = 2 * radius ^ 2;

switch lower(char(mode))
    case '3d'
        vertexCoords = vertices;
        sourceCoords = activeXYZ;
    case 'yz'
        vertexCoords = vertices(:, 2:3);
        sourceCoords = activeXYZ(:, 2:3);
end

sourceNorm = sum(sourceCoords .^ 2, 2).';
for firstIdx = 1:chunkSize:size(vertices, 1)
    lastIdx = min(firstIdx + chunkSize - 1, size(vertices, 1));
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
