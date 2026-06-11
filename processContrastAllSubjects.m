%% processContrastAllSubjects.m
% Run subject-level sLORETA contrasts for every subject.
%
% This is a template batch script. Before running it, set dataFile and
% dataVariable below to match your own full group data file.
%
% The contrast is always:
%
%   condition 1 - condition 2
%
% Arrange your input EEG matrix so condition 1 is the positive condition and
% condition 2 is the negative condition.

clear;
clc;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);

dataFile = 'yourGroupData.mat';
dataVariable = 'yourDataVariable';
condition1Idx = 1;
condition2Idx = 2;
timeIndex = 317;

if strcmp(dataFile, 'yourGroupData.mat') || strcmp(dataVariable, 'yourDataVariable')
    error('Edit dataFile and dataVariable near the top of this script before running the batch example.');
end

outputDir = fullfile('outputs', sprintf('contrast_cond%02d_minus_cond%02d_time%03d', ...
    condition1Idx, condition2Idx, timeIndex));

loaded = load(dataFile, dataVariable);
allData = loaded.(dataVariable);
nSubjects = size(allData, 4);

fprintf('Found %d subjects in %s.\n', nSubjects, dataFile);
fprintf('Contrast: condition %d - condition %d.\n', condition1Idx, condition2Idx);
fprintf('Using time point %d for every subject.\n', timeIndex);

allOutputs = repmat(struct( ...
    'imageFile', '', ...
    'csvFile', '', ...
    'matFile', '', ...
    'timeIndex', [], ...
    'sourceContrast', [], ...
    'sourceXYZ', [], ...
    'regionContrast', [], ...
    'metadata', []), nSubjects, 1);

for subjectIdx = 1:nSubjects
    outputPrefix = sprintf('rewp_sub%02d_cond%02dminus%02d_t%03d', ...
        subjectIdx, condition1Idx, condition2Idx, timeIndex);

    fprintf('\nSubject %d of %d: %s\n', subjectIdx, nSubjects, outputPrefix);

    allOutputs(subjectIdx) = doContrastLoreta(dataFile, ...
        'DataVariable', dataVariable, ...
        'SubjectIdx', subjectIdx, ...
        'Condition1Idx', condition1Idx, ...
        'Condition2Idx', condition2Idx, ...
        'TimeIndex', timeIndex, ...
        'OutputDir', outputDir, ...
        'OutputPrefix', outputPrefix);
end

summaryFile = fullfile(outputDir, sprintf('contrast_cond%02d_minus_cond%02d_time%03d_summary.mat', ...
    condition1Idx, condition2Idx, timeIndex));
save(summaryFile, 'allOutputs', 'condition1Idx', 'condition2Idx', ...
    'timeIndex', 'nSubjects', '-v7');

fprintf('\nDone processing contrast condition %d - condition %d.\n', ...
    condition1Idx, condition2Idx);
fprintf('Subject-level contrast files are in: %s\n', outputDir);
fprintf('Batch summary saved: %s\n', summaryFile);
