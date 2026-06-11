%% exampleUseage.m
% This file shows the simplest way to run sLORETA plotting for one subject
% and one condition.
%
% You can run this whole file by pressing the Run button in MATLAB.
%
% The main function is:
%
%   doPlotsLoreta(...)
%
% It does three things:
%
%   1. Runs sLORETA on one subject and one condition.
%   2. Saves a six-view brain image as a PNG file.
%   3. Saves a CSV file that ranks brain regions from most active to least
%      active.
%
% It also saves a .mat file that can be used later for group-level analysis.

clear;
clc;

%% Step 1: Go to the folder that contains these sLORETA files
% This lets the example work no matter where the toolbox folder is located.

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);

%% Step 2: Load the input data
% The current beginner example data file is:
%
%   exampleSubjectData.mat
%
% It contains a variable called:
%
%   exampleData
%
% The expected data shape is:
%
%   channels x time x conditions
%
% For the current file:
%
%   63 channels x 400 time points x 2 conditions
%
% Data loading is deliberately done outside doPlotsLoreta. That keeps the
% main function independent of your file names and variable names.

load('exampleSubjectData.mat', 'exampleData');

%% Step 3: Choose which subject and condition to analyze
% MATLAB uses 1-based indexing.
%
% This example file already contains one extracted subject, so subjectIdx is
% mainly saved as metadata in the output files.
% Condition 1 means the first condition in the file.

subjectIdx = 1;
conditionIdx = 1;

%% Step 4: Choose the time point to plot
% Leave this empty [] if you want the function to automatically pick the
% time point with the strongest average source activity.
%
% Or set a specific time point, for example:
%
%   timeIndex = 317;
%
% For group-level analysis, it is usually better to use the same time point
% for every subject.

timeIndex = [];

%% Step 5: Choose where the output files should go
% The folder will be created if it does not already exist.

outputDir = 'outputs';

%% Step 6: Choose a name for the output files
% The function will add endings such as:
%
%   _sixview.png
%   _regions.csv
%   _source.mat
%
% For example, if outputPrefix is 'rewp_sub01_cond01', the image file will
% be:
%
%   rewp_sub01_cond01_sixview.png

outputPrefix = 'rewp_sub01_cond01';

%% Step 7: Run sLORETA and make the outputs
% Most users only need to edit the variables above this line.
%
% MRI slice notes:
%
%   By default, the MRI slice plot makes a 2 x 3 montage. The first row
%   shows sagittal/coronal/axial slices through the strongest activity above
%   z = -35 mm. The second row shows a second nearby-but-distinct peak.
%   This avoids very low slices where the template MRI can show face/neck
%   anatomy instead of useful brain anatomy.
%
%   To force a specific MRI/source slice, add:
%
%       'SliceXYZ', [0 0 30], ...
%
%   To allow the raw maximum activity location even if it is very low, add:
%
%       'AutoSliceMinZ', -Inf, ...
%
%   To make the MRI activity overlay smoother or sharper, adjust:
%
%       'OverlaySigmaMm', 10, ...
%
%   To include a little more or less natural MRI edge around the brain mask,
%   adjust:
%
%       'BrainMaskDilateVoxels', 2, ...
%
%   To show only one MRI row, add:
%
%       'NumSliceRows', 1, ...
%
%   Figures stay open by default. For batch processing, add:
%
%       'CloseFigures', true, ...
%
% Inputs explained:
%
%   exampleData
%       The already-loaded EEG/ERP matrix.
%
%   'SubjectIdx'
%       Which subject to analyze. For this example file, there is only one
%       subject already extracted, so this value is just saved as metadata.
%
%   'ConditionIdx'
%       Which condition to analyze.
%
%   'TimeIndex'
%       Which time point to plot. Use [] for automatic selection.
%
%   'OutputDir'
%       Folder where the PNG, CSV, and MAT files are saved.
%
%   'OutputPrefix'
%       Beginning of each output filename.

outputs = doPlotsLoreta(exampleData, ...
    'SubjectIdx', subjectIdx, ...
    'ConditionIdx', conditionIdx, ...
    'TimeIndex', timeIndex, ...
    'OutputDir', outputDir, ...
    'OutputPrefix', outputPrefix);

%% Step 8: Show the saved filenames in the MATLAB Command Window

fprintf('\nDone.\n');
fprintf('Brain image saved here:\n%s\n\n', outputs.imageFile);
fprintf('Internal slice image saved here:\n%s\n\n', outputs.sliceFile);
fprintf('MRI slice image saved here:\n%s\n\n', outputs.mriSliceFile);
fprintf('Ranked brain-region CSV saved here:\n%s\n\n', outputs.csvFile);
fprintf('Group-analysis MAT file saved here:\n%s\n\n', outputs.matFile);
fprintf('Time point plotted: %d\n', outputs.timeIndex);

%% Step 9: Optional notes
% The CSV file ranks regions using the DKT atlas by default.
%
% The .mat file contains:
%
%   sourceERP
%       The source-localized data. This is useful for later group analysis.
%
%   sourceValues
%       The source values at the plotted time point.
%
%   regionRanking
%       The same region-ranking information saved in the CSV file.
%
%   metadata
%       Information about the subject, condition, time point, and settings.
%
% To analyze a different subject or condition, change subjectIdx or
% conditionIdx near the top of this file and run the script again.
