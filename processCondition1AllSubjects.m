%% processCondition1AllSubjects.m
% Run sLORETA for every subject in condition 1.
%
% This is a template batch script. Before running it, set dataFile and
% dataVariable below to match your own full group data file.
%
% This script creates one output set per subject:
%
%   *_sixview.png     brain image
%   *_regions.csv    ranked brain regions
%   *_source.mat     source data for group-level analysis
%
% The group-level analysis should use the *_source.mat files.

clear;
clc;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);

dataFile = 'yourGroupData.mat';
dataVariable = 'yourDataVariable';
conditionIdx = 1;

if strcmp(dataFile, 'yourGroupData.mat') || strcmp(dataVariable, 'yourDataVariable')
    error('Edit dataFile and dataVariable near the top of this script before running the batch example.');
end

% Choose the time point to use for every subject.
%
% Use a number for group-consistent analysis, for example:
%
%   timeIndex = 317;
%
% Use [] only if you want each subject to use their own strongest time point.
% For group-level analysis, a fixed time point is usually easier to interpret.
timeIndex = 317;

if isempty(timeIndex)
    outputDir = fullfile('outputs', 'condition01_autoPeak');
    timeLabel = 'auto';
else
    outputDir = fullfile('outputs', sprintf('condition01_time%03d', timeIndex));
    timeLabel = sprintf('t%03d', timeIndex);
end

loaded = load(dataFile, dataVariable);
allData = loaded.(dataVariable);
nSubjects = size(allData, 4);

fprintf('Found %d subjects in %s.\n', nSubjects, dataFile);
fprintf('Processing condition %d.\n', conditionIdx);
if isempty(timeIndex)
    fprintf('Using each subject''s strongest time point.\n');
else
    fprintf('Using time point %d for every subject.\n', timeIndex);
end

allOutputs = repmat(struct( ...
    'imageFile', '', ...
    'csvFile', '', ...
    'matFile', '', ...
    'timeIndex', [], ...
    'sourceERP', [], ...
    'sourceXYZ', [], ...
    'regionRanking', [], ...
    'metadata', []), nSubjects, 1);

for subjectIdx = 1:nSubjects
    outputPrefix = sprintf('rewp_sub%02d_cond%02d_%s', ...
        subjectIdx, conditionIdx, timeLabel);

    fprintf('\nSubject %d of %d: %s\n', subjectIdx, nSubjects, outputPrefix);

    allOutputs(subjectIdx) = doPlotsLoreta(dataFile, ...
        'DataVariable', dataVariable, ...
        'SubjectIdx', subjectIdx, ...
        'ConditionIdx', conditionIdx, ...
        'TimeIndex', timeIndex, ...
        'OutputDir', outputDir, ...
        'OutputPrefix', outputPrefix);
end

summaryFile = fullfile(outputDir, sprintf('condition%02d_batch_summary.mat', conditionIdx));
save(summaryFile, 'allOutputs', 'conditionIdx', 'nSubjects', 'timeIndex', '-v7');

fprintf('\nDone processing condition %d.\n', conditionIdx);
fprintf('Subject-level files are in: %s\n', outputDir);
fprintf('Batch summary saved: %s\n', summaryFile);
