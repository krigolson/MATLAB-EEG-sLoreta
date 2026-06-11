function outputs = doPlotsLoreta(dataInput, varargin)
%DOPLOTSLORETA Run subject-level sLORETA and save plot/CSV/MAT outputs.
%
%   outputs = doPlotsLoreta(dataInput)
%   outputs = doPlotsLoreta(dataInput, 'Name', value, ...)
%
%   dataInput can be:
%       - numeric channels x time, channels x time x conditions, or
%         channels x time x conditions x subjects
%       - a .mat filename containing the data variable
%
%   Common example:
%
%       outputs = doPlotsLoreta('exampleSubjectData.mat', ...
%           'DataVariable', 'exampleData', ...
%           'SubjectIdx', 1, ...
%           'ConditionIdx', 1, ...
%           'OutputDir', 'outputs');
%
%   Saved outputs:
%       - six-view cortex PNG
%       - ranked brain-region CSV
%       - .mat file with sourceERP/sourceXYZ/metadata for group analysis

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
addParameter(parser, 'ConditionIdx', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(parser, 'TimeIndex', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1));
addParameter(parser, 'Lambda', 0.05, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(parser, 'Reference', 'average', @(x) ischar(x) || isstring(x));
addParameter(parser, 'AtlasName', 'DKT', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ThresholdPercentile', 84, @(x) isnumeric(x) && ...
    isscalar(x) && x >= 0 && x <= 100);
addParameter(parser, 'ColorPercentile', 99, @(x) isnumeric(x) && ...
    isscalar(x) && x >= 0 && x <= 100);
addParameter(parser, 'ProjectionRadius', 35, @(x) isnumeric(x) && ...
    isscalar(x) && isfinite(x) && x > 0);
addParameter(parser, 'ProjectionMode', '3d', @(x) ischar(x) || isstring(x));
addParameter(parser, 'MaxSources', Inf, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(parser, 'Region', 'all', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ACCBox', [-15 15 0 50 0 55], @(x) isnumeric(x) && numel(x) == 6);
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

[data, dataInfo] = load_or_select_data(dataInput, opts);

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
validate_channel_match(data, model, channelLabels, leadfieldPath);

requestedDataSize = [size(data, 1), size(data, 2), 1, 1];
data = reshape(data, requestedDataSize);

inverse = sloreta_make_inverse(model.leadfield, ...
    'SourceXYZ', model.sourceXYZ, ...
    'Lambda', opts.Lambda, ...
    'Reference', char(opts.Reference), ...
    'Orientation', 'free', ...
    'Standardize', true);

sourceERP = sloreta_apply(inverse, data, 'Output', 'norm');
sourceERP = reshape(sourceERP, [inverse.nSources, size(data, 2), 1, 1]);

source2d = squeeze(sourceERP);
timeIndex = opts.TimeIndex;
if isempty(timeIndex)
    [~, timeIndex] = max(mean(abs(source2d), 1));
else
    timeIndex = round(timeIndex);
end

if isempty(char(opts.OutputPrefix))
    outputPrefix = sprintf('sLoreta_sub%02d_cond%02d_t%03d', ...
        opts.SubjectIdx, opts.ConditionIdx, timeIndex);
else
    outputPrefix = char(opts.OutputPrefix);
end

matFile = fullfile(outputDir, [outputPrefix '_source.mat']);
imageFile = fullfile(outputDir, [outputPrefix '_sixview.png']);
sliceFile = fullfile(outputDir, [outputPrefix '_slices.png']);
mriSliceFile = fullfile(outputDir, [outputPrefix '_mri_slices.png']);
csvFile = fullfile(outputDir, [outputPrefix '_regions.csv']);

sourceValues = abs(source2d(:, timeIndex));
regionRanking = rank_cortex_regions(sourceValues, inverse.sourceXYZ, opts);
write_region_csv(csvFile, regionRanking);

plotSourceFile = matFile;
metadata = struct();
metadata.createdBy = mfilename;
metadata.dataInfo = dataInfo;
metadata.subjectIdx = opts.SubjectIdx;
metadata.conditionIdx = opts.ConditionIdx;
metadata.timeIndex = timeIndex;
metadata.lambda = opts.Lambda;
metadata.reference = char(opts.Reference);
metadata.atlasName = char(opts.AtlasName);
metadata.thresholdPercentile = opts.ThresholdPercentile;
metadata.colorPercentile = opts.ColorPercentile;
metadata.projectionRadius = opts.ProjectionRadius;
metadata.projectionMode = char(opts.ProjectionMode);
metadata.region = char(opts.Region);
metadata.autoSliceMinZ = opts.AutoSliceMinZ;
metadata.overlaySigmaMm = opts.OverlaySigmaMm;
metadata.numSliceRows = opts.NumSliceRows;
metadata.minPeakDistanceMm = opts.MinPeakDistanceMm;
metadata.brainMaskLabel = opts.BrainMaskLabel;
metadata.brainMaskDilateVoxels = opts.BrainMaskDilateVoxels;
metadata.outputPrefix = outputPrefix;
metadata.requestedDataSize = requestedDataSize;
metadata.sourceSize = [inverse.nSources, size(data, 2), 1, 1];

save(matFile, 'sourceERP', 'sourceValues', 'inverse', 'channelLabels', ...
    'metadata', 'regionRanking', '-v7');

plot_sloreta_cortex(plotSourceFile, ...
    'TimeIndex', timeIndex, ...
    'View', 'six', ...
    'ThresholdPercentile', opts.ThresholdPercentile, ...
    'ColorPercentile', opts.ColorPercentile, ...
    'ProjectionRadius', opts.ProjectionRadius, ...
    'ProjectionMode', char(opts.ProjectionMode), ...
    'MaxSources', opts.MaxSources, ...
    'Region', char(opts.Region), ...
    'ACCBox', opts.ACCBox, ...
    'OutputPng', imageFile);

if logical(opts.MakeSlicePlot)
    sliceFig = doPlotLoretaSlices(matFile, ...
        'MapType', 'activity', ...
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
        'MapType', 'activity', ...
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
outputs.sourceERP = sourceERP;
outputs.sourceXYZ = inverse.sourceXYZ;
outputs.regionRanking = regionRanking;
outputs.metadata = metadata;

fprintf('sLORETA plot saved: %s\n', imageFile);
if logical(opts.MakeSlicePlot)
    fprintf('sLORETA slice plot saved: %s\n', sliceFile);
end
if logical(opts.MakeMriSlicePlot)
    fprintf('sLORETA MRI slice plot saved: %s\n', mriSliceFile);
end
fprintf('Region CSV saved: %s\n', csvFile);
fprintf('Group-analysis MAT saved: %s\n', matFile);

end

function [data, info] = load_or_select_data(dataInput, opts)
info = struct();
if isnumeric(dataInput)
    raw = dataInput;
    info.source = 'numeric input';
else
    dataFile = char(dataInput);
    loaded = load(dataFile, char(opts.DataVariable));
    varName = char(opts.DataVariable);
    if ~isfield(loaded, varName)
        error([mfilename ':MissingDataVariable'], ...
            'Could not find variable %s in %s.', varName, dataFile);
    end
    raw = loaded.(varName);
    info.source = dataFile;
    info.variable = varName;
end

if ndims(raw) == 2
    data = raw;
elseif ndims(raw) == 3
    data = raw(:, :, opts.ConditionIdx);
else
    data = raw(:, :, opts.ConditionIdx, opts.SubjectIdx);
end

info.rawSize = size(raw);
info.selectedSize = size(data);
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

function ranking = rank_cortex_regions(sourceValues, sourceXYZ, opts)
cortexFile = fullfile(fileparts(mfilename('fullpath')), ...
    'templates', 'cortex', 'brainstorm_icbm152_cortex_pial_low.mat');
cortex = load(cortexFile, 'Vertices', 'Faces', 'Atlas');
vertices = cortex.Vertices;
if max(abs(vertices(:))) < 1
    vertices = vertices * 1000;
end

sourceValues = abs(sourceValues(:));
mask = source_region_mask(sourceXYZ, char(opts.Region), opts.ACCBox);
sourceValues(~mask) = 0;
sourceValues = keep_strongest_sources(sourceValues, opts.MaxSources);

vertexValues = project_active_sources_to_vertices(vertices, sourceXYZ, sourceValues, ...
    opts.ProjectionRadius, char(opts.ProjectionMode));
vertexValues = smooth_vertex_values(cortex.Faces, vertexValues, 3);

atlasIdx = find(strcmpi({cortex.Atlas.Name}, char(opts.AtlasName)), 1);
if isempty(atlasIdx)
    available = strjoin({cortex.Atlas.Name}, ', ');
    error([mfilename ':BadAtlas'], ...
        'AtlasName %s not found. Available atlases: %s', ...
        char(opts.AtlasName), available);
end

scouts = cortex.Atlas(atlasIdx).Scouts;
nScouts = numel(scouts);
rows = repmat(struct('rank', [], 'atlas', '', 'label', '', 'region', '', ...
    'meanActivity', [], 'maxActivity', [], 'sumActivity', [], 'nVertices', []), ...
    nScouts, 1);

for idx = 1:nScouts
    scoutVertices = scouts(idx).Vertices(:);
    scoutVertices = scoutVertices(scoutVertices >= 1 & scoutVertices <= numel(vertexValues));
    vals = vertexValues(scoutVertices);
    rows(idx).atlas = cortex.Atlas(atlasIdx).Name;
    rows(idx).label = scouts(idx).Label;
    rows(idx).region = scouts(idx).Region;
    rows(idx).meanActivity = mean(vals);
    rows(idx).maxActivity = max(vals);
    rows(idx).sumActivity = sum(vals);
    rows(idx).nVertices = numel(vals);
end

[~, order] = sort([rows.meanActivity], 'descend');
rows = rows(order);
for idx = 1:numel(rows)
    rows(idx).rank = idx;
end

ranking = rows;
end

function write_region_csv(csvFile, ranking)
ensure_parent_dir(csvFile);
fid = fopen(csvFile, 'w');
if fid < 0
    error([mfilename ':CannotWriteCSV'], 'Could not write %s.', csvFile);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'rank,atlas,label,region,meanActivity,maxActivity,sumActivity,nVertices\n');
for idx = 1:numel(ranking)
    fprintf(fid, '%d,%s,%s,%s,%.12g,%.12g,%.12g,%d\n', ...
        ranking(idx).rank, ...
        csv_escape(ranking(idx).atlas), ...
        csv_escape(ranking(idx).label), ...
        csv_escape(ranking(idx).region), ...
        ranking(idx).meanActivity, ...
        ranking(idx).maxActivity, ...
        ranking(idx).sumActivity, ...
        ranking(idx).nVertices);
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

function text = csv_escape(value)
text = char(value);
text = strrep(text, '"', '""');
text = ['"', text, '"'];
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
        error([mfilename ':BadRegion'], ...
            'Region must be ''all'' or ''acc''.');
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

switch lower(char(mode))
    case '3d'
        vertexCoords = vertices;
        sourceCoords = activeXYZ;
    case 'yz'
        vertexCoords = vertices(:, 2:3);
        sourceCoords = activeXYZ(:, 2:3);
    otherwise
        error([mfilename ':BadProjectionMode'], ...
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
