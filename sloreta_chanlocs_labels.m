function labels = sloreta_chanlocs_labels(chanlocs)
%SLORETA_CHANLOCS_LABELS Return channel labels from an EEGLAB-like struct.
%
%   labels = sloreta_chanlocs_labels(chanlocs)
%
%   This helper has no EEGLAB dependency. It only expects a struct array with
%   a labels field, like the chanlocs saved in matlocs.mat.

if ~isstruct(chanlocs) || ~isfield(chanlocs, 'labels')
    error('%s:BadChanlocs', mfilename, ...
        'chanlocs must be a struct array with a labels field.');
end

labels = cell(numel(chanlocs), 1);
for idx = 1:numel(chanlocs)
    labels{idx} = chanlocs(idx).labels;
end

end
