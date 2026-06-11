function P = sloreta_average_reference_projector(nChannels)
%SLORETA_AVERAGE_REFERENCE_PROJECTOR Project channel data to average reference.
%
%   P = sloreta_average_reference_projector(nChannels)
%
%   Returns an nChannels-by-nChannels matrix that removes the across-channel
%   mean at each time point. Apply the same projector to both scalp data and
%   the lead field.

validateattributes(nChannels, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, mfilename, 'nChannels');

P = eye(nChannels) - ones(nChannels) ./ nChannels;

end
