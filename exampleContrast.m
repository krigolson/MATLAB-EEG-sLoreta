%% exampleContrast.m
% Run a subject-level sLORETA contrast for one subject.
%
% The contrast is always:
%
%   condition 1 - condition 2
%
% The labels do not matter to the code. Condition 1 might be wins, targets,
% congruent trials, etc. Condition 2 might be losses, standards, incongruent
% trials, etc. Arrange your EEG matrix so condition 1 is the positive
% condition.

clear;
clc;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);

% Use the bundled sample data.
%
% sampleGrandERP is:
%
%   channels x time x conditions x subjects
%
% The contrast function will compare condition1Idx against condition2Idx.
% The code does not care what those conditions mean.
dataFile = 'sampleGrandERP.mat';
dataVariable = 'sampleGrandERP';
chanlocsFile = 'sampleGrandERP.mat';

subjectIdx = 1;
condition1Idx = 1;
condition2Idx = 2;
timeIndex = 317;

outputDir = fullfile('outputs', 'contrast_single_subject');
outputPrefix = sprintf('rewp_sub%02d_cond%02dminus%02d_t%03d', ...
    subjectIdx, condition1Idx, condition2Idx, timeIndex);

outputs = doContrastLoreta(dataFile, ...
    'DataVariable', dataVariable, ...
    'ChanlocsFile', chanlocsFile, ...
    'SubjectIdx', subjectIdx, ...
    'Condition1Idx', condition1Idx, ...
    'Condition2Idx', condition2Idx, ...
    'TimeIndex', timeIndex, ...
    'OutputDir', outputDir, ...
    'OutputPrefix', outputPrefix);

fprintf('\nDone.\n');
fprintf('Signed contrast image:\n%s\n\n', outputs.imageFile);
fprintf('Signed contrast slice image:\n%s\n\n', outputs.sliceFile);
fprintf('Signed contrast MRI slice image:\n%s\n\n', outputs.mriSliceFile);
fprintf('Signed contrast region CSV:\n%s\n\n', outputs.csvFile);
fprintf('Signed contrast MAT file:\n%s\n\n', outputs.matFile);
