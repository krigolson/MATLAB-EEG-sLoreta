function sourceData = sloreta_apply(inverse, data, varargin)
%SLORETA_APPLY Apply a prepared sLORETA inverse to channel data.
%
%   sourceData = sloreta_apply(inverse, data)
%   sourceData = sloreta_apply(inverse, data, 'Output', 'norm')
%
%   data can be channels x time x conditions, or channels x any trailing
%   dimensions. For the current REwP file this means both subject-level
%   channels x time x conditions and group arrays like
%   channels x time x conditions x subjects are supported.
%
%   Output options:
%       'components'   dipoles x trailing dimensions. Default.
%       'norm'         sources x trailing dimensions for free-orientation
%                      models; same as components for scalar models.

parser = inputParser;
parser.FunctionName = mfilename;
addRequired(parser, 'inverse', @isstruct);
addRequired(parser, 'data', @(x) validateattributes(x, {'numeric'}, ...
    {'real', 'finite', 'nonempty'}, mfilename, 'data'));
addParameter(parser, 'Output', 'components', @(x) ischar(x) || isstring(x));
parse(parser, inverse, data, varargin{:});

output = lower(char(parser.Results.Output));

requiredFields = {'W', 'projector', 'nChannels', 'nSources', 'nOrient'};
for idx = 1:numel(requiredFields)
    if ~isfield(inverse, requiredFields{idx})
        error('%s:BadInverse', mfilename, ...
            'Inverse struct is missing field ''%s''.', requiredFields{idx});
    end
end

dataSize = size(data);
if dataSize(1) ~= inverse.nChannels
    error('%s:ChannelMismatch', mfilename, ...
        'Data has %d channels, but inverse expects %d.', ...
        dataSize(1), inverse.nChannels);
end

trailingSize = dataSize(2:end);
data2d = reshape(double(data), inverse.nChannels, []);
data2d = inverse.projector * data2d;

source2d = inverse.W * data2d;

switch output
    case 'components'
        sourceData = reshape(source2d, [size(source2d, 1), trailingSize]);
    case 'norm'
        if inverse.nOrient == 1
            sourceData = reshape(source2d, [inverse.nSources, trailingSize]);
        else
            source3d = reshape(source2d, inverse.nOrient, inverse.nSources, []);
            source2d = squeeze(sqrt(sum(source3d .^ 2, 1)));
            sourceData = reshape(source2d, [inverse.nSources, trailingSize]);
        end
    otherwise
        error('%s:BadOutput', mfilename, ...
            'Output must be ''components'' or ''norm''.');
end

end
