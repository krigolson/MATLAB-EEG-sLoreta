function inverse = sloreta_make_inverse(leadfield, varargin)
%SLORETA_MAKE_INVERSE Build a standardized minimum-norm inverse operator.
%
%   inverse = sloreta_make_inverse(leadfield)
%   inverse = sloreta_make_inverse(leadfield, 'Name', value, ...)
%
%   Required input
%   --------------
%   leadfield : channels x dipoles numeric matrix. For fixed/scalar source
%       orientation, dipoles equals sources. For free orientation, dipoles
%       should be 3*sources with x/y/z columns grouped by source.
%
%   Name-value options
%   ------------------
%   'Lambda'           Relative Tikhonov regularization. Default: 0.05.
%                      The absolute ridge is Lambda * trace(L*L')/channels.
%   'Reference'        'average' or 'none'. Default: 'average'.
%   'Orientation'      'scalar', 'free', or 'auto'. Default: 'auto'.
%   'SourceXYZ'        Optional sources x 3 coordinates.
%   'Standardize'      true/false. Default: true.
%
%   Output
%   ------
%   inverse is a struct with fields:
%       W              dipoles x channels inverse operator
%       leadfield      referenced leadfield used to build W
%       projector      channel reference projector
%       lambda         relative lambda
%       ridge          absolute ridge
%       nChannels      number of channels
%       nDipoles       number of inverse rows
%       nSources       number of source locations
%       nOrient        1 or 3
%       sourceXYZ      optional coordinates
%
%   Notes
%   -----
%   This function intentionally does not depend on EEGLAB. If chanlocs come
%   from EEGLAB, use them only to keep channel order aligned with leadfield.

parser = inputParser;
parser.FunctionName = mfilename;
addRequired(parser, 'leadfield', @(x) validateattributes(x, {'numeric'}, ...
    {'2d', 'real', 'finite', 'nonempty'}, mfilename, 'leadfield'));
addParameter(parser, 'Lambda', 0.05, @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, 'Lambda'));
addParameter(parser, 'Reference', 'average', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Orientation', 'auto', @(x) ischar(x) || isstring(x));
addParameter(parser, 'SourceXYZ', [], @(x) isempty(x) || isnumeric(x));
addParameter(parser, 'Standardize', true, @(x) islogical(x) && isscalar(x));
parse(parser, leadfield, varargin{:});

L = double(parser.Results.leadfield);
lambda = parser.Results.Lambda;
reference = lower(char(parser.Results.Reference));
orientation = lower(char(parser.Results.Orientation));
sourceXYZ = parser.Results.SourceXYZ;
doStandardize = parser.Results.Standardize;

[nChannels, nDipoles] = size(L);

switch reference
    case 'average'
        projector = sloreta_average_reference_projector(nChannels);
    case 'none'
        projector = eye(nChannels);
    otherwise
        error('%s:BadReference', mfilename, ...
            'Reference must be ''average'' or ''none''.');
end

L = projector * L;

if isempty(sourceXYZ)
    if strcmp(orientation, 'free')
        if mod(nDipoles, 3) ~= 0
            error('%s:BadOrientation', mfilename, ...
                'Free-orientation leadfields must have 3*sources columns.');
        end
        nOrient = 3;
        nSources = nDipoles / 3;
    else
        nOrient = 1;
        nSources = nDipoles;
    end
else
    validateattributes(sourceXYZ, {'numeric'}, ...
        {'2d', 'real', 'finite', 'ncols', 3}, mfilename, 'SourceXYZ');
    nSources = size(sourceXYZ, 1);
    if nDipoles == nSources
        nOrient = 1;
    elseif nDipoles == 3 * nSources
        nOrient = 3;
    else
        error('%s:BadSourceXYZ', mfilename, ...
            'Leadfield columns must equal sources or 3*sources.');
    end
end

if strcmp(orientation, 'scalar') && nOrient ~= 1
    error('%s:BadOrientation', mfilename, ...
        'Orientation is scalar, but leadfield/source dimensions imply free orientation.');
elseif strcmp(orientation, 'free') && nOrient ~= 3
    error('%s:BadOrientation', mfilename, ...
        'Orientation is free, but leadfield/source dimensions do not imply 3 orientations.');
elseif ~ismember(orientation, {'auto', 'scalar', 'free'})
    error('%s:BadOrientation', mfilename, ...
        'Orientation must be ''auto'', ''scalar'', or ''free''.');
end

gram = L * L.';
scale = trace(gram) / max(nChannels, 1);
ridge = lambda * scale;
systemMatrix = gram + ridge * eye(nChannels);

% Minimum-norm inverse. Right matrix division avoids explicitly inverting
% the channel covariance-like system.
W = L.' / systemMatrix;

if doStandardize
    W = standardize_inverse(W, L, nSources, nOrient);
end

inverse = struct();
inverse.W = W;
inverse.leadfield = L;
inverse.projector = projector;
inverse.lambda = lambda;
inverse.ridge = ridge;
inverse.nChannels = nChannels;
inverse.nDipoles = nDipoles;
inverse.nSources = nSources;
inverse.nOrient = nOrient;
inverse.sourceXYZ = sourceXYZ;
inverse.reference = reference;
inverse.standardized = doStandardize;

end

function Wstd = standardize_inverse(W, L, nSources, nOrient)
resolution = W * L;
Wstd = W;

if nOrient == 1
    denom = sqrt(max(real(diag(resolution)), eps));
    Wstd = W ./ denom;
    return;
end

for sourceIdx = 1:nSources
    rows = (sourceIdx - 1) * 3 + (1:3);
    block = real((resolution(rows, rows) + resolution(rows, rows).') ./ 2);
    [vectors, values] = eig(block);
    values = max(diag(values), eps);
    invSqrt = vectors * diag(1 ./ sqrt(values)) * vectors.';
    Wstd(rows, :) = invSqrt * W(rows, :);
end

end
