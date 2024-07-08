function [ data ] = fir_filterdcpadded(b, a, data, causal, usefftfilt)

% Defaults
if nargin < 5 || isempty(usefftfilt)
    usefftfilt = 0;
end
if nargin < 4 || isempty(causal)
    causal = 0;
end

% Check arguments
if nargin < 3
    error('Not enough input arguments.');
end

% Is FIR?
if ~isscalar(a) || a ~= 1
    error('Not a FIR filter. onepass-zerophase and onepass-minphase filtering is available for FIR filters only.')
end

% Group delay
if mod(length(b), 2) ~= 1
    error('Filter order is not even.');
end
groupDelay = (length(b) - 1) / 2;

% Filter symmetry
isSym = all(b(1:groupDelay) == b(end:-1:groupDelay + 2));
isAntisym = all([b(1:groupDelay) == -b(end:-1:groupDelay + 2) b(groupDelay + 1) == 0]);
if causal == 0 && ~(isSym || isAntisym)
    error('Filter is not anti-/symmetric. For onepass-zerophase filtering the filter must be anti-/symmetric.')
end

% Padding
if causal
    startPad = repmat(data(1, :), [2 * groupDelay 1]);
    endPad = [];
else
    startPad = repmat(data(1, :), [groupDelay 1]);
    endPad = repmat(data(end, :), [groupDelay 1]);
end

% Filter data (with double precision)
isSingle = isa(data, 'single');

if usefftfilt
    data = fftfilt(double(b), double([startPad; data; endPad]));
else
    data = filter(double(b), 1, double([startPad; data; endPad])); % Pad and filter with double precision
end

% Convert to single
if isSingle
    data = single(data);
end

% Remove padded data
data = data(2 * groupDelay + 1:end, :);

end