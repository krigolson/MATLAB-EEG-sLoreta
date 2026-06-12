function outputs = doGroupContrastLoreta(inputDir, varargin)
%DOGROUPCONTRASTLORETA Group test for condition1-condition2 sLORETA contrasts.
%
%   outputs = doGroupContrastLoreta(inputDir)
%
%   inputDir should contain *_contrast.mat files created by doContrastLoreta.
%   The group test is a one-sample t-test of signed regional contrast values
%   against zero.

if nargin < 1 || isempty(inputDir)
    inputDir = fullfile('outputs', 'contrast_cond01_minus_cond02_time317');
end

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'FilePattern', '*_contrast.mat', @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputDir', inputDir, @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputPrefix', 'group_contrast_cond01_minus_cond02_t317', @(x) ischar(x) || isstring(x));
addParameter(parser, 'AtlasName', 'DKT', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Measure', 'meanContrast', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(parser, 'Correction', 'fdr', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Tail', 'both', @(x) ischar(x) || isstring(x));
addParameter(parser, 'MapStatistic', 't', @(x) ischar(x) || isstring(x));
addParameter(parser, 'BrainTemplate', 'brainnet', @(x) ischar(x) || isstring(x));
addParameter(parser, 'SurfaceFile', '', @(x) ischar(x) || isstring(x));
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
inputDir = char(inputDir);
outputDir = char(opts.OutputDir);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

files = dir(fullfile(inputDir, char(opts.FilePattern)));
if isempty(files)
    error([mfilename ':NoFiles'], ...
        'No files matched %s in %s.', char(opts.FilePattern), inputDir);
end

[regionLabels, regionCodes, atlasName, values, subjectFiles] = ...
    load_contrast_matrix(inputDir, files, char(opts.Measure));
[sourceValues, inverse] = load_source_matrix(subjectFiles);

stats = compute_region_stats(values, char(opts.Tail));
stats.pCorrected = correct_pvalues(stats.p, char(opts.Correction));
stats.significant = stats.pCorrected <= opts.Alpha;
groupSourceStats = compute_source_stats(sourceValues, char(opts.Tail), ...
    char(opts.Correction), opts.Alpha);

regionStats = make_region_stats_struct(regionLabels, regionCodes, atlasName, values, stats);
significantStats = regionStats([regionStats.significant]);
significantStats = sort_region_stats(significantStats, char(opts.MapStatistic));

csvFile = fullfile(outputDir, [char(opts.OutputPrefix) '_all_regions.csv']);
matFile = fullfile(outputDir, [char(opts.OutputPrefix) '_group_contrast.mat']);
imageFile = fullfile(outputDir, [char(opts.OutputPrefix) '_significant_regions.png']);
pMapFile = fullfile(outputDir, [char(opts.OutputPrefix) '_uncorrected_pmap.png']);
sliceFile = fullfile(outputDir, [char(opts.OutputPrefix) '_group_slices.png']);
mriSliceFile = fullfile(outputDir, [char(opts.OutputPrefix) '_group_mri_slices.png']);

write_group_contrast_csv(csvFile, regionStats);
plot_group_contrast_regions(significantStats, atlasName, char(opts.MapStatistic), ...
    imageFile, logical(opts.CloseFigures), opts.BrainTemplate, opts.SurfaceFile);
plot_group_contrast_regions(regionStats, atlasName, 'signedLogP', ...
    pMapFile, logical(opts.CloseFigures), opts.BrainTemplate, opts.SurfaceFile);

metadata = struct();
metadata.createdBy = mfilename;
metadata.inputDir = inputDir;
metadata.filePattern = char(opts.FilePattern);
metadata.subjectFiles = subjectFiles;
metadata.nSubjects = numel(subjectFiles);
metadata.atlasName = atlasName;
metadata.measure = char(opts.Measure);
metadata.alpha = opts.Alpha;
metadata.correction = char(opts.Correction);
metadata.tail = char(opts.Tail);
metadata.mapStatistic = char(opts.MapStatistic);
metadata.contrast = 'condition1_minus_condition2';
metadata.autoSliceMinZ = opts.AutoSliceMinZ;
metadata.overlaySigmaMm = opts.OverlaySigmaMm;
metadata.numSliceRows = opts.NumSliceRows;
metadata.minPeakDistanceMm = opts.MinPeakDistanceMm;
metadata.brainMaskLabel = opts.BrainMaskLabel;
metadata.brainMaskDilateVoxels = opts.BrainMaskDilateVoxels;

save(matFile, 'regionStats', 'significantStats', 'values', 'regionLabels', ...
    'regionCodes', 'sourceValues', 'groupSourceStats', 'inverse', ...
    'metadata', '-v7');

if logical(opts.MakeSlicePlot)
    sliceFig = doPlotLoretaSlices(matFile, ...
        'MapType', 'meanContrast', ...
        'SliceXYZ', opts.SliceXYZ, ...
        'SliceX', opts.SliceX, ...
        'SliceY', opts.SliceY, ...
        'SliceZ', opts.SliceZ, ...
        'OutputPng', sliceFile);
    if logical(opts.CloseFigures)
        close(sliceFig.figure);
    end
end

if logical(opts.MakeMriSlicePlot)
    mriSliceFig = doPlotLoretaMriSlices(matFile, ...
        'MapType', 'meanContrast', ...
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
        'OutputPng', mriSliceFile);
    if logical(opts.CloseFigures)
        close(mriSliceFig.figure);
    end
end

outputs = struct();
outputs.imageFile = imageFile;
outputs.pMapFile = pMapFile;
outputs.sliceFile = sliceFile;
outputs.mriSliceFile = mriSliceFile;
outputs.csvFile = csvFile;
outputs.matFile = matFile;
outputs.regionStats = regionStats;
outputs.significantStats = significantStats;
outputs.groupSourceStats = groupSourceStats;
outputs.metadata = metadata;

fprintf('Group contrast significant-region map saved: %s\n', imageFile);
fprintf('Group contrast uncorrected p-map saved: %s\n', pMapFile);
if logical(opts.MakeSlicePlot)
    fprintf('Group contrast slice plot saved: %s\n', sliceFile);
end
if logical(opts.MakeMriSlicePlot)
    fprintf('Group contrast MRI slice plot saved: %s\n', mriSliceFile);
end
fprintf('Group contrast all-region CSV saved: %s\n', csvFile);
fprintf('Group contrast MAT saved: %s\n', matFile);
fprintf('Significant contrast regions: %d of %d\n', numel(significantStats), numel(regionStats));

end

function [sourceValues, inverse] = load_source_matrix(subjectFiles)
nSubjects = numel(subjectFiles);
sourceValues = [];
inverse = [];

for subjIdx = 1:nSubjects
    loaded = load(subjectFiles{subjIdx}, 'sourceValues', 'inverse');
    if ~isfield(loaded, 'sourceValues') || ~isfield(loaded, 'inverse')
        error([mfilename ':MissingSourceValues'], ...
            '%s must contain sourceValues and inverse.', subjectFiles{subjIdx});
    end
    if subjIdx == 1
        nSources = numel(loaded.sourceValues);
        sourceValues = zeros(nSubjects, nSources);
        inverse = loaded.inverse;
    elseif numel(loaded.sourceValues) ~= size(sourceValues, 2)
        error([mfilename ':SourceSizeMismatch'], ...
            'Source count differs in %s.', subjectFiles{subjIdx});
    end
    sourceValues(subjIdx, :) = loaded.sourceValues(:).';
end
end

function groupSourceStats = compute_source_stats(sourceValues, tail, correction, alpha)
stats = compute_region_stats(sourceValues, tail);
pCorrected = correct_pvalues(stats.p, correction);

groupSourceStats = struct();
groupSourceStats.n = stats.n;
groupSourceStats.df = stats.df;
groupSourceStats.mean = stats.mean(:);
groupSourceStats.sd = stats.sd(:);
groupSourceStats.sem = stats.sem(:);
groupSourceStats.t = stats.t(:);
groupSourceStats.p = stats.p(:);
groupSourceStats.pCorrected = pCorrected(:);
groupSourceStats.cohenDz = stats.cohenDz(:);
groupSourceStats.significant = groupSourceStats.pCorrected <= alpha;
groupSourceStats.correction = correction;
groupSourceStats.tail = tail;
groupSourceStats.alpha = alpha;
end

function [labels, regionCodes, atlasName, values, subjectFiles] = load_contrast_matrix(inputDir, files, measure)
nSubjects = numel(files);
labels = {};
regionCodes = {};
atlasName = '';
values = [];
subjectFiles = cell(nSubjects, 1);

for subjIdx = 1:nSubjects
    subjectFiles{subjIdx} = fullfile(inputDir, files(subjIdx).name);
    loaded = load(subjectFiles{subjIdx}, 'regionContrast');
    rr = loaded.regionContrast;
    if subjIdx == 1
        labels = {rr.label}.';
        regionCodes = {rr.region}.';
        atlasName = rr(1).atlas;
        values = zeros(nSubjects, numel(rr));
    end

    currentLabels = {rr.label}.';
    for regionIdx = 1:numel(labels)
        matchIdx = find(strcmp(labels{regionIdx}, currentLabels), 1);
        if isempty(matchIdx)
            error([mfilename ':RegionMismatch'], ...
                'Could not find region %s in %s.', labels{regionIdx}, files(subjIdx).name);
        end
        values(subjIdx, regionIdx) = rr(matchIdx).(measure);
    end
end
end

function stats = compute_region_stats(values, tail)
n = size(values, 1);
df = n - 1;
meanValue = mean(values, 1);
sdValue = std(values, 0, 1);
semValue = sdValue ./ sqrt(n);
tValue = meanValue ./ semValue;
tValue(sdValue == 0 & meanValue > 0) = Inf;
tValue(sdValue == 0 & meanValue == 0) = 0;
tValue(sdValue == 0 & meanValue < 0) = -Inf;

stats = struct();
stats.n = n;
stats.df = df;
stats.mean = meanValue;
stats.sd = sdValue;
stats.sem = semValue;
stats.t = tValue;
stats.p = t_pvalue(tValue, df, tail);
stats.cohenDz = meanValue ./ sdValue;
stats.cohenDz(sdValue == 0) = Inf;
end

function p = t_pvalue(t, df, tail)
x = df ./ (df + t .^ 2);
twoTail = betainc(x, df / 2, 0.5);
rightTail = 0.5 .* twoTail;
rightTail(t < 0) = 1 - 0.5 .* twoTail(t < 0);
switch lower(tail)
    case 'right'
        p = rightTail;
    case 'left'
        p = 1 - rightTail;
    case 'both'
        p = twoTail;
    otherwise
        error([mfilename ':BadTail'], 'Tail must be right, left, or both.');
end
p(isnan(p)) = 1;
end

function pCorrected = correct_pvalues(p, correction)
switch lower(correction)
    case 'none'
        pCorrected = p;
    case 'fdr'
        pCorrected = fdr_bh(p);
    otherwise
        error([mfilename ':BadCorrection'], 'Correction must be fdr or none.');
end
end

function adjusted = fdr_bh(p)
p = p(:);
m = numel(p);
[sortedP, order] = sort(p, 'ascend');
adjustedSorted = sortedP .* m ./ (1:m).';
for idx = m-1:-1:1
    adjustedSorted(idx) = min(adjustedSorted(idx), adjustedSorted(idx + 1));
end
adjustedSorted = min(adjustedSorted, 1);
adjusted = zeros(size(p));
adjusted(order) = adjustedSorted;
adjusted = adjusted.';
end

function rows = make_region_stats_struct(labels, regionCodes, atlasName, values, stats)
nRegions = numel(labels);
rows = repmat(struct('rank', [], 'atlas', '', 'label', '', 'region', '', ...
    'nSubjects', [], 'meanContrast', [], 'sdContrast', [], 'semContrast', [], ...
    't', [], 'df', [], 'p', [], 'pCorrected', [], 'cohenDz', [], ...
    'significant', [], 'direction', ''), nRegions, 1);

for idx = 1:nRegions
    rows(idx).rank = idx;
    rows(idx).atlas = atlasName;
    rows(idx).label = labels{idx};
    rows(idx).region = regionCodes{idx};
    rows(idx).nSubjects = size(values, 1);
    rows(idx).meanContrast = stats.mean(idx);
    rows(idx).sdContrast = stats.sd(idx);
    rows(idx).semContrast = stats.sem(idx);
    rows(idx).t = stats.t(idx);
    rows(idx).df = stats.df;
    rows(idx).p = stats.p(idx);
    rows(idx).pCorrected = stats.pCorrected(idx);
    rows(idx).cohenDz = stats.cohenDz(idx);
    rows(idx).significant = stats.significant(idx);
    rows(idx).direction = contrast_direction(stats.mean(idx));
end

rows = sort_region_stats(rows, 'absT');
for idx = 1:numel(rows)
    rows(idx).rank = idx;
end
end

function sortedStats = sort_region_stats(regionStats, mapStatistic)
if isempty(regionStats)
    sortedStats = regionStats;
    return;
end
switch lower(mapStatistic)
    case 't'
        sortValues = [regionStats.t];
    case 'abst'
        sortValues = abs([regionStats.t]);
    case {'mean', 'meancontrast'}
        sortValues = [regionStats.meanContrast];
    case {'absmean', 'absmeancontrast'}
        sortValues = abs([regionStats.meanContrast]);
    otherwise
        error([mfilename ':BadMapStatistic'], ...
        'MapStatistic must be t, absT, meanContrast, or absMeanContrast.');
end
[~, order] = sort(sortValues, 'descend');
sortedStats = regionStats(order);
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

function write_group_contrast_csv(csvFile, rows)
fid = fopen(csvFile, 'w');
if fid < 0
    error([mfilename ':CannotWriteCSV'], 'Could not write %s.', csvFile);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, ['rank,atlas,label,region,nSubjects,meanContrast,sdContrast,', ...
    'semContrast,t,df,p,pCorrected,cohenDz,significant,direction\n']);
for idx = 1:numel(rows)
    fprintf(fid, '%d,%s,%s,%s,%d,%.12g,%.12g,%.12g,%.12g,%d,%.12g,%.12g,%.12g,%d,%s\n', ...
        idx, csv_escape(rows(idx).atlas), csv_escape(rows(idx).label), ...
        csv_escape(rows(idx).region), rows(idx).nSubjects, rows(idx).meanContrast, ...
        rows(idx).sdContrast, rows(idx).semContrast, rows(idx).t, rows(idx).df, ...
        rows(idx).p, rows(idx).pCorrected, rows(idx).cohenDz, ...
        rows(idx).significant, csv_escape(rows(idx).direction));
end
end

function plot_group_contrast_regions(rows, atlasName, mapStatistic, outputPng, closeFigures, brainTemplate, surfaceFile)
cortexFile = fullfile(fileparts(mfilename('fullpath')), ...
    'templates', 'cortex', 'brainstorm_icbm152_cortex_pial_low.mat');
cortex = load(cortexFile, 'Vertices', 'Faces', 'Atlas');
vertices = cortex.Vertices;
faces = cortex.Faces;
if max(abs(vertices(:))) < 1
    vertices = vertices * 1000;
end
vertexValues = zeros(size(vertices, 1), 1);
atlasIdx = find(strcmpi({cortex.Atlas.Name}, atlasName), 1);

for idx = 1:numel(rows)
    scoutIdx = find(strcmp({cortex.Atlas(atlasIdx).Scouts.Label}, rows(idx).label), 1);
    if isempty(scoutIdx)
        continue;
    end
    scoutVertices = cortex.Atlas(atlasIdx).Scouts(scoutIdx).Vertices(:);
    scoutVertices = scoutVertices(scoutVertices >= 1 & scoutVertices <= numel(vertexValues));
    vertexValues(scoutVertices) = region_map_value(rows(idx), mapStatistic);
end

fig = plot_sloreta_signed_cortex(vertices, vertexValues, ...
    'View', 'six', ...
    'ProjectionRadius', 1, ...
    'BrainTemplate', brainTemplate, ...
    'SurfaceFile', surfaceFile, ...
    'OutputPng', outputPng);
if closeFigures
    close(fig);
end
end

function value = region_map_value(row, mapStatistic)
switch lower(mapStatistic)
    case {'t', 'abst'}
        value = row.t;
    case {'mean', 'meancontrast', 'absmean', 'absmeancontrast'}
        value = row.meanContrast;
    case {'signedlogp', 'pmap'}
        value = sign(row.meanContrast) * -log10(max(row.p, realmin));
    otherwise
        error([mfilename ':BadMapStatistic'], ...
            'MapStatistic must be t, absT, meanContrast, absMeanContrast, or signedLogP.');
end
end

function text = csv_escape(value)
text = char(value);
text = strrep(text, '"', '""');
text = ['"', text, '"'];
end
