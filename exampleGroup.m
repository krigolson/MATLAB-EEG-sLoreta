%% exampleGroup.m
% Run a group-level sLORETA contrast analysis.
%
% This is the recommended group workflow.
%
% It assumes you have already run:
%
%   processContrastAllSubjects
%
% That script creates one *_contrast.mat file per subject.
%
% The group test asks:
%
%   Is condition 1 - condition 2 reliably different from zero?
%
% Positive values mean:
%
%   condition 1 > condition 2
%
% Negative values mean:
%
%   condition 2 > condition 1

clear;
clc;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);

inputDir = fullfile('outputs', 'contrast_cond01_minus_cond02_time317');
outputDir = fullfile('outputs', 'group_contrast_cond01_minus_cond02_time317');
outputPrefix = 'group_cond01_minus_cond02_t317';

outputs = doGroupContrastLoreta(inputDir, ...
    'OutputDir', outputDir, ...
    'OutputPrefix', outputPrefix, ...
    'Alpha', 0.05, ...
    'Correction', 'fdr', ...
    'Tail', 'both', ...
    'AtlasName', 'DKT', ...
    'Measure', 'meanContrast', ...
    'MapStatistic', 't');

fprintf('\nDone.\n');
fprintf('Significant contrast-region brain image:\n%s\n\n', outputs.imageFile);
fprintf('Uncorrected contrast p-map image:\n%s\n\n', outputs.pMapFile);
fprintf('Group contrast slice image:\n%s\n\n', outputs.sliceFile);
fprintf('Group contrast MRI slice image:\n%s\n\n', outputs.mriSliceFile);
fprintf('All-region contrast statistics CSV:\n%s\n\n', outputs.csvFile);
fprintf('Group contrast MAT file:\n%s\n\n', outputs.matFile);
fprintf('Number of significant contrast regions: %d\n', numel(outputs.significantStats));
