function leadfieldFile = make_leadfield_from_dipfit_template(varargin)
%MAKE_LEADFIELD_FROM_DIPFIT_TEMPLATE Build leadfield.mat from local templates.
%
%   make_leadfield_from_dipfit_template
%   make_leadfield_from_dipfit_template('Resolution', 10)
%
%   This uses FieldTrip to compute a template EEG lead field from local
%   copies of EEGLAB DIPFIT's standard BEM head model and standard 10-05
%   electrode file. The output leadfield.mat is independent of EEGLAB and
%   FieldTrip and can be used by sloreta_make_inverse.
%
%   The saved variables are:
%       leadfield       channels x 3*sources
%       sourceXYZ       sources x 3, in mm
%       channelLabels   channels x 1 cell array
%       leadfieldInfo   provenance/settings struct

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'Resolution', 10, @(x) validateattributes(x, ...
    {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'Resolution'));
addParameter(parser, 'FieldTripPath', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ChanlocsFile', 'matlocs.mat', @(x) ischar(x) || isstring(x));
addParameter(parser, 'DipfitPath', ...
    fullfile(fileparts(mfilename('fullpath')), 'templates', 'dipfit'), ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'HeadmodelFile', 'standard_vol.mat', ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'ElectrodeFile', 'standard_1005.elc', ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputFile', 'leadfield.mat', ...
    @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});

resolution = parser.Results.Resolution;
fieldtripPath = char(parser.Results.FieldTripPath);
chanlocsFile = char(parser.Results.ChanlocsFile);
dipfitPath = char(parser.Results.DipfitPath);
headmodelName = char(parser.Results.HeadmodelFile);
electrodeName = char(parser.Results.ElectrodeFile);
outputFile = char(parser.Results.OutputFile);

restorePath = path;
cleanup = onCleanup(@() path(restorePath));

rootDir = fileparts(mfilename('fullpath'));
if ~is_absolute_path(chanlocsFile)
    chanlocsFile = fullfile(rootDir, chanlocsFile);
end
if ~is_absolute_path(outputFile)
    outputFile = fullfile(rootDir, outputFile);
end

if ~isempty(fieldtripPath)
    addpath(fieldtripPath);
end
if exist('ft_defaults', 'file') ~= 2
    error([mfilename ':MissingFieldTrip'], ...
        ['FieldTrip is needed only to regenerate leadfield.mat. Add ', ...
        'FieldTrip to the MATLAB path or call this function with ', ...
        '''FieldTripPath'', ''/path/to/fieldtrip''.']);
end
ft_defaults;

load(chanlocsFile, 'chanlocs');
channelLabels = sloreta_chanlocs_labels(chanlocs);

headmodelFile = fullfile(dipfitPath, 'standard_BEM', headmodelName);
elecFile = fullfile(dipfitPath, 'standard_BEM', 'elec', electrodeName);

if ~exist(headmodelFile, 'file')
    error([mfilename ':MissingHeadmodel'], ...
        'Could not find head model file: %s', headmodelFile);
end
if ~exist(elecFile, 'file')
    error([mfilename ':MissingElectrodes'], ...
        'Could not find electrode file: %s', elecFile);
end

fprintf('Loading head model: %s\n', headmodelFile);
headmodelStruct = load(headmodelFile);
headmodel = ft_convert_units(headmodelStruct.vol, 'mm');

fprintf('Loading electrodes: %s\n', elecFile);
elec = ft_read_sens(elecFile, 'senstype', 'eeg');
elec = ft_convert_units(elec, 'mm');
elec = select_and_order_electrodes(elec, channelLabels);

cfg = [];
cfg.headmodel = headmodel;
cfg.elec = elec;
cfg.channel = channelLabels;
cfg.resolution = resolution;
cfg.unit = 'mm';
cfg.reducerank = 'no';
cfg.feedback = 'text';

fprintf('Preparing lead field at %.2f mm source-grid resolution...\n', resolution);
sourcemodel = ft_prepare_leadfield(cfg);

inside = sourcemodel.inside;
if islogical(inside)
    insideIdx = find(inside(:));
elseif isnumeric(inside) && all(inside(:) == 0 | inside(:) == 1)
    insideIdx = find(inside(:));
elseif isnumeric(inside)
    insideIdx = inside(:);
else
    error([mfilename ':BadInsideField'], ...
        'Unexpected sourcemodel.inside type.');
end

nChannels = numel(channelLabels);
nSources = numel(insideIdx);
leadfield = zeros(nChannels, 3 * nSources);
sourceXYZ = sourcemodel.pos(insideIdx, :);

for sourceIdx = 1:nSources
    lf = sourcemodel.leadfield{insideIdx(sourceIdx)};
    if isempty(lf)
        error([mfilename ':EmptyLeadfield'], ...
            'Empty lead field at source index %d.', insideIdx(sourceIdx));
    end
    if size(lf, 1) ~= nChannels || size(lf, 2) ~= 3
        error([mfilename ':BadLeadfieldShape'], ...
            'Expected %d x 3 lead field, got %d x %d.', ...
            nChannels, size(lf, 1), size(lf, 2));
    end
    columns = (sourceIdx - 1) * 3 + (1:3);
    leadfield(:, columns) = lf;
end

leadfieldInfo = struct();
leadfieldInfo.source = 'Local copy of EEGLAB DIPFIT standard_BEM + FieldTrip ft_prepare_leadfield';
leadfieldInfo.headmodelFile = headmodelFile;
leadfieldInfo.elecFile = elecFile;
leadfieldInfo.fieldtripPath = fieldtripPath;
leadfieldInfo.chanlocsFile = chanlocsFile;
leadfieldInfo.resolutionMm = resolution;
leadfieldInfo.nChannels = nChannels;
leadfieldInfo.nSources = nSources;
leadfieldInfo.orientation = 'free';
leadfieldInfo.createdBy = mfilename;

save(outputFile, 'leadfield', 'sourceXYZ', 'channelLabels', 'leadfieldInfo', '-v7.3');

fprintf('Saved %s with leadfield size [%s] and sourceXYZ size [%s].\n', ...
    outputFile, num2str(size(leadfield)), num2str(size(sourceXYZ)));

leadfieldFile = outputFile;

end

function tf = is_absolute_path(pathIn)
pathIn = char(pathIn);
tf = startsWith(pathIn, filesep) || ~isempty(regexp(pathIn, '^[A-Za-z]:[\\/]', 'once'));
end

function elecOut = select_and_order_electrodes(elecIn, channelLabels)
elecLabels = elecIn.label(:);
nChannels = numel(channelLabels);
selection = zeros(nChannels, 1);

for idx = 1:nChannels
    match = find(strcmpi(channelLabels{idx}, elecLabels), 1);
    if isempty(match)
        error([mfilename ':MissingElectrode'], ...
            'Could not find channel %s in the standard electrode file.', ...
            channelLabels{idx});
    end
    selection(idx) = match;
end

elecOut = elecIn;
elecOut.label = elecIn.label(selection);

if isfield(elecIn, 'chanpos')
    elecOut.chanpos = elecIn.chanpos(selection, :);
end
if isfield(elecIn, 'elecpos')
    elecOut.elecpos = elecIn.elecpos(selection, :);
end
if isfield(elecIn, 'tra')
    elecOut.tra = eye(nChannels);
end

end
