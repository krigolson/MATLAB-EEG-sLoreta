function outputs = doMakeLeadfield(varargin)
%DOMAKELEADFIELD Build a montage-specific leadfield file.
%
%   outputs = doMakeLeadfield
%   outputs = doMakeLeadfield('ChanlocsFile', 'matlocs_32chan.mat')
%   outputs = doMakeLeadfield('ChanlocsFile', 'matlocs_32chan.mat', ...
%       'OutputFile', 'leadfield_32chan.mat', ...
%       'FieldTripPath', '/path/to/fieldtrip')
%
%   The channel labels and order in ChanlocsFile define the rows of the
%   saved leadfield. Your EEG matrix must use the same channel order.

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'ChanlocsFile', 'sampleGrandERP.mat', @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputFile', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'FieldTripPath', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Resolution', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(parser, 'Force', false, @(x) islogical(x) || isnumeric(x));
parse(parser, varargin{:});
opts = parser.Results;

rootDir = fileparts(mfilename('fullpath'));
chanlocsFile = resolve_path(char(opts.ChanlocsFile), rootDir);
loaded = load(chanlocsFile, 'chanlocs');
if ~isfield(loaded, 'chanlocs')
    error([mfilename ':MissingChanlocs'], ...
        'ChanlocsFile must contain a variable named chanlocs.');
end

channelLabels = sloreta_chanlocs_labels(loaded.chanlocs);
nChannels = numel(channelLabels);

outputFile = char(opts.OutputFile);
if isempty(outputFile)
    if ismember(lower(get_filename(chanlocsFile)), {'matlocs.mat', 'samplegranderp.mat'}) && nChannels == 63
        outputFile = 'leadfield.mat';
    else
        outputFile = sprintf('leadfield_%dchan.mat', nChannels);
    end
end
outputFile = resolve_path(outputFile, rootDir);

if exist(outputFile, 'file') && ~logical(opts.Force)
    fprintf('Leadfield already exists: %s\n', outputFile);
    fprintf('Use ''Force'', true to rebuild it.\n');
else
    make_leadfield_from_dipfit_template( ...
        'ChanlocsFile', chanlocsFile, ...
        'OutputFile', outputFile, ...
        'FieldTripPath', char(opts.FieldTripPath), ...
        'Resolution', opts.Resolution);
end

outputs = struct();
outputs.leadfieldFile = outputFile;
outputs.chanlocsFile = chanlocsFile;
outputs.nChannels = nChannels;
outputs.channelLabels = channelLabels;
outputs.resolution = opts.Resolution;

fprintf('Leadfield file ready: %s\n', outputFile);
fprintf('Channels: %d\n', nChannels);

end

function pathOut = resolve_path(pathIn, rootDir)
pathIn = char(pathIn);
if is_absolute_path(pathIn)
    pathOut = pathIn;
else
    pathOut = fullfile(rootDir, pathIn);
end
end

function tf = is_absolute_path(pathIn)
pathIn = char(pathIn);
tf = startsWith(pathIn, filesep) || ~isempty(regexp(pathIn, '^[A-Za-z]:[\\/]', 'once'));
end

function name = get_filename(pathIn)
[~, base, ext] = fileparts(pathIn);
name = [base ext];
end
