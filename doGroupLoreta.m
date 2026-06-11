function outputs = doGroupLoreta(inputDir, varargin)
%DOGROUPLORETA Group-level statistics and significant-region brain map.
%
%   outputs = doGroupLoreta(inputDir)
%   outputs = doGroupLoreta(inputDir, 'Name', value, ...)
%
%   inputDir should contain subject-level *_source.mat files created by
%   doPlotsLoreta. The group test is run across subjects for each atlas
%   region. Only statistically significant regions are shown on the output
%   brain map and written to the CSV file.
%
%   Important note:
%   The current subject-level values are sLORETA norm/absolute activity
%   values. A one-sample test against zero is useful as a first pipeline
%   check, but a future condition-difference or baseline-corrected analysis
%   will usually be more scientifically meaningful.

if nargin < 1 || isempty(inputDir)
    inputDir = fullfile('outputs', 'condition01_time317');
end

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'FilePattern', '*_source.mat', @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputDir', inputDir, @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputPrefix', 'group_condition01_time317', @(x) ischar(x) || isstring(x));
addParameter(parser, 'AtlasName', 'DKT', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Measure', 'meanActivity', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(parser, 'Correction', 'fdr', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Tail', 'right', @(x) ischar(x) || isstring(x));
addParameter(parser, 'TestValue', 0, @(x) isnumeric(x) && isscalar(x));
addParameter(parser, 'MapStatistic', 't', @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});

opts = parser.Results;
inputDir = char(inputDir);
outputDir = char(opts.OutputDir);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

files = dir(fullfile(inputDir, char(opts.FilePattern)));
if isempty(files)
    error('%s:NoFiles', mfilename, ...
        'No files matched %s in %s.', char(opts.FilePattern), inputDir);
end

[regionLabels, regionCodes, atlasName, values, subjectFiles] = ...
    load_region_matrix(inputDir, files, char(opts.Measure), char(opts.AtlasName));

stats = compute_region_stats(values, opts.TestValue, char(opts.Tail));
stats.pCorrected = correct_pvalues(stats.p, char(opts.Correction));
stats.significant = stats.pCorrected <= opts.Alpha;

regionStats = make_region_stats_struct(regionLabels, regionCodes, atlasName, ...
    values, stats, opts);
significantStats = regionStats([regionStats.significant]);
significantStats = sort_region_stats(significantStats, char(opts.MapStatistic));

if numel(significantStats) == numel(regionStats)
    fprintf(['Note: all regions are significant. This often happens when ', ...
        'testing nonnegative sLORETA norm values against zero. Consider ', ...
        'using condition differences or baseline-corrected values for a ', ...
        'more selective group map.\n']);
end

csvFile = fullfile(outputDir, [char(opts.OutputPrefix) '_all_regions.csv']);
matFile = fullfile(outputDir, [char(opts.OutputPrefix) '_group.mat']);
imageFile = fullfile(outputDir, [char(opts.OutputPrefix) '_significant_regions.png']);
pMapFile = fullfile(outputDir, [char(opts.OutputPrefix) '_uncorrected_pmap.png']);

write_group_csv(csvFile, regionStats);
plot_group_significant_regions(significantStats, char(opts.AtlasName), ...
    char(opts.MapStatistic), imageFile);
plot_group_significant_regions(regionStats, char(opts.AtlasName), ...
    'negLogP', pMapFile);

metadata = struct();
metadata.createdBy = mfilename;
metadata.inputDir = inputDir;
metadata.filePattern = char(opts.FilePattern);
metadata.subjectFiles = subjectFiles;
metadata.nSubjects = numel(subjectFiles);
metadata.atlasName = char(opts.AtlasName);
metadata.measure = char(opts.Measure);
metadata.alpha = opts.Alpha;
metadata.correction = char(opts.Correction);
metadata.tail = char(opts.Tail);
metadata.testValue = opts.TestValue;
metadata.mapStatistic = char(opts.MapStatistic);

save(matFile, 'regionStats', 'significantStats', 'values', 'regionLabels', ...
    'regionCodes', 'metadata', '-v7');

outputs = struct();
outputs.imageFile = imageFile;
outputs.pMapFile = pMapFile;
outputs.csvFile = csvFile;
outputs.matFile = matFile;
outputs.regionStats = regionStats;
outputs.significantStats = significantStats;
outputs.metadata = metadata;

fprintf('Group significant-region map saved: %s\n', imageFile);
fprintf('Group uncorrected p-map saved: %s\n', pMapFile);
fprintf('Group all-region CSV saved: %s\n', csvFile);
fprintf('Group MAT file saved: %s\n', matFile);
fprintf('Significant regions: %d of %d\n', numel(significantStats), numel(regionStats));

end

function [labels, regionCodes, atlasName, values, subjectFiles] = ...
    load_region_matrix(inputDir, files, measure, atlasRequested)
nSubjects = numel(files);
labels = {};
regionCodes = {};
atlasName = atlasRequested;
values = [];
subjectFiles = cell(nSubjects, 1);

for subjIdx = 1:nSubjects
    subjectFiles{subjIdx} = fullfile(inputDir, files(subjIdx).name);
    loaded = load(subjectFiles{subjIdx}, 'regionRanking');
    rr = loaded.regionRanking;
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
            error('%s:RegionMismatch', mfilename, ...
                'Could not find region %s in %s.', labels{regionIdx}, files(subjIdx).name);
        end
        if ~isfield(rr(matchIdx), measure)
            error('%s:BadMeasure', mfilename, ...
                'Measure %s not found in regionRanking.', measure);
        end
        values(subjIdx, regionIdx) = rr(matchIdx).(measure);
    end
end

if ~strcmpi(atlasName, atlasRequested)
    warning('%s:AtlasMismatch', mfilename, ...
        'Requested atlas %s, but files report atlas %s.', atlasRequested, atlasName);
end
end

function stats = compute_region_stats(values, testValue, tail)
n = size(values, 1);
df = n - 1;
diffValues = values - testValue;
meanValue = mean(values, 1);
sdValue = std(values, 0, 1);
semValue = sdValue ./ sqrt(n);
tValue = mean(diffValues, 1) ./ semValue;
tValue(sdValue == 0 & mean(diffValues, 1) > 0) = Inf;
tValue(sdValue == 0 & mean(diffValues, 1) == 0) = 0;
tValue(sdValue == 0 & mean(diffValues, 1) < 0) = -Inf;

pValue = t_pvalue(tValue, df, tail);
cohenDz = mean(diffValues, 1) ./ sdValue;
cohenDz(sdValue == 0) = Inf;

stats = struct();
stats.n = n;
stats.df = df;
stats.mean = meanValue;
stats.sd = sdValue;
stats.sem = semValue;
stats.t = tValue;
stats.p = pValue;
stats.cohenDz = cohenDz;
end

function p = t_pvalue(t, df, tail)
t = double(t);
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
        error('%s:BadTail', mfilename, 'Tail must be right, left, or both.');
end

p(isnan(p)) = 1;
p(t == Inf & strcmpi(tail, 'right')) = 0;
p(t == -Inf & strcmpi(tail, 'left')) = 0;
end

function pCorrected = correct_pvalues(p, correction)
switch lower(correction)
    case 'none'
        pCorrected = p;
    case 'fdr'
        pCorrected = fdr_bh(p);
    otherwise
        error('%s:BadCorrection', mfilename, ...
            'Correction must be fdr or none.');
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

function regionStats = make_region_stats_struct(labels, regionCodes, atlasName, values, stats, opts)
nRegions = numel(labels);
regionStats = repmat(struct('rank', [], 'atlas', '', 'label', '', 'region', '', ...
    'nSubjects', [], 'meanActivity', [], 'sdActivity', [], 'semActivity', [], ...
    't', [], 'df', [], 'p', [], 'pCorrected', [], 'cohenDz', [], ...
    'significant', []), nRegions, 1);

for idx = 1:nRegions
    regionStats(idx).rank = idx;
    regionStats(idx).atlas = atlasName;
    regionStats(idx).label = labels{idx};
    regionStats(idx).region = regionCodes{idx};
    regionStats(idx).nSubjects = size(values, 1);
    regionStats(idx).meanActivity = stats.mean(idx);
    regionStats(idx).sdActivity = stats.sd(idx);
    regionStats(idx).semActivity = stats.sem(idx);
    regionStats(idx).t = stats.t(idx);
    regionStats(idx).df = stats.df;
    regionStats(idx).p = stats.p(idx);
    regionStats(idx).pCorrected = stats.pCorrected(idx);
    regionStats(idx).cohenDz = stats.cohenDz(idx);
    regionStats(idx).significant = stats.significant(idx);
end

regionStats = sort_region_stats(regionStats, char(opts.MapStatistic));
for idx = 1:numel(regionStats)
    regionStats(idx).rank = idx;
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
    case {'mean', 'meanactivity'}
        sortValues = [regionStats.meanActivity];
    case {'cohendz', 'effect'}
        sortValues = [regionStats.cohenDz];
    otherwise
        error('%s:BadMapStatistic', mfilename, ...
            'MapStatistic must be t, meanActivity, or cohenDz.');
end

[~, order] = sort(sortValues, 'descend');
sortedStats = regionStats(order);
end

function write_group_csv(csvFile, statsRows)
fid = fopen(csvFile, 'w');
if fid < 0
    error('%s:CannotWriteCSV', mfilename, 'Could not write %s.', csvFile);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, ['rank,atlas,label,region,nSubjects,meanActivity,sdActivity,', ...
    'semActivity,t,df,p,pCorrected,cohenDz,significant\n']);
for idx = 1:numel(statsRows)
    fprintf(fid, '%d,%s,%s,%s,%d,%.12g,%.12g,%.12g,%.12g,%d,%.12g,%.12g,%.12g,%d\n', ...
        idx, ...
        csv_escape(statsRows(idx).atlas), ...
        csv_escape(statsRows(idx).label), ...
        csv_escape(statsRows(idx).region), ...
        statsRows(idx).nSubjects, ...
        statsRows(idx).meanActivity, ...
        statsRows(idx).sdActivity, ...
        statsRows(idx).semActivity, ...
        statsRows(idx).t, ...
        statsRows(idx).df, ...
        statsRows(idx).p, ...
        statsRows(idx).pCorrected, ...
        statsRows(idx).cohenDz, ...
        statsRows(idx).significant);
end
end

function text = csv_escape(value)
text = char(value);
text = strrep(text, '"', '""');
text = ['"', text, '"'];
end

function plot_group_significant_regions(significantStats, atlasName, mapStatistic, outputPng)
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
if isempty(atlasIdx)
    error('%s:BadAtlas', mfilename, 'Atlas %s not found.', atlasName);
end

for idx = 1:numel(significantStats)
    scoutIdx = find(strcmp({cortex.Atlas(atlasIdx).Scouts.Label}, significantStats(idx).label), 1);
    if isempty(scoutIdx)
        continue;
    end
    verticesInScout = cortex.Atlas(atlasIdx).Scouts(scoutIdx).Vertices(:);
    verticesInScout = verticesInScout(verticesInScout >= 1 & verticesInScout <= numel(vertexValues));
    vertexValues(verticesInScout) = max(vertexValues(verticesInScout), ...
        region_map_value(significantStats(idx), mapStatistic));
end

colors = significant_region_colors(vertexValues);
fig = figure('Color', 'k', 'Name', 'Group significant sLORETA regions');
set(fig, 'Renderer', 'opengl', 'InvertHardcopy', 'off', 'Position', [80 80 1500 900]);

viewNames = {'left', 'top', 'right', 'front', 'back', 'bottom'};
for viewIdx = 1:numel(viewNames)
    ax = subplot(2, 3, viewIdx, 'Parent', fig);
    render_cortex_axes(ax, vertices, faces, colors, viewNames{viewIdx}, 1.02);
    title(ax, viewNames{viewIdx}, 'Color', 'w', 'FontWeight', 'normal');
end

save_figure_png(fig, outputPng, 250);
end

function value = region_map_value(row, mapStatistic)
switch lower(mapStatistic)
    case 't'
        value = row.t;
    case {'mean', 'meanactivity'}
        value = row.meanActivity;
    case {'cohendz', 'effect'}
        value = row.cohenDz;
    case {'neglogp', 'pmap'}
        value = -log10(max(row.p, realmin));
end
end

function colors = significant_region_colors(values)
baseGray = repmat([0.58 0.58 0.56], numel(values), 1);
active = values > 0;
colors = baseGray;
if ~any(active)
    return;
end

activeValues = values(active);
lo = min(activeValues);
hi = max(activeValues);
if hi <= lo
    hi = lo + eps;
end
scaled = (activeValues - lo) ./ (hi - lo);
scaled = min(max(scaled, 0), 1);
hotMap = hot(256);
idx = max(1, min(256, round(1 + scaled * 255)));
colors(active, :) = hotMap(idx, :);
end

function render_cortex_axes(ax, vertices, faces, colors, viewSpec, zoomFactor)
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
view(ax, view_to_az_el(viewSpec));
camlight(ax, 'headlight');
camlight(ax, -80, 20);
camlight(ax, 100, 25);
lighting(ax, 'gouraud');
material(ax, 'dull');
camzoom(ax, zoomFactor);
end

function azel = view_to_az_el(viewSpec)
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
