function [index, d] = photobleach_index(signal, sigma, threshold)
% Finds the point in a FRET donor/acceptor signal where 
% photobleaching occurs.
%
% Inputs
% ------
% signal : (Tx1)
% 	A 1D donor or acceptor signal 
% sigma : float
% 	Optional width (std dev) of gaussian kernel for smoothing of 
%   detection signal (default = 1.0)
% threshold : float (optional)
% 	Return 0 if no peak found above a certain threshold
% 
% Outputs
% ------- 
% index : int
%	Index of the bleaching point. Returns length(signal) if no 
%   bleaching point is found.
% d : (Tx1)
%	Detection signal 
% 
% Jan-Willem van de Meent
% $Revision: 1.30 $  $Date: 2012/02/13$
% $Revision: 1.20 $  $Date: 2011/05/04$
% $Revision: 1.10 $  $Date: 2011/04/27$
% $Revision: 1.00 $  $Date: 2011/04/15$

T = length(signal);

% (profiled out)
% for t = 1:length(signal)
%     % calculate mean and std dev
%     % of signal from start to t 
%     mf(t) = mean(signal(1:t));
%     sf(t) = std(signal(1:t));
%     % calculate mean and std dev
%     % of signal from t to end
%     mb(t) = mean(signal(t:end));
%     sb(t) = std(signal(t:end));
% end

% forward mean: mf(t) = mean(signal(1:t))
mf = (tril(ones(T)) * signal) ./ (1:T)';
% foward std dev: std(signal(1:t))
sf = sqrt((tril(ones(T)) * signal.^2) ./ (1:T)' - mf.^2);
% backward mean: mb(t) = mean(signal(t:end))
mb = (triu(ones(T)) * signal) ./ (T:-1:1)';
%backward std dev: std(signal(t:end))
sb = sqrt((triu(ones(T)) * signal.^2) ./ (T:-1:1)' - mb.^2);

% create gaussian for time smoothing
if nargin < 2
    sigma = 5;
end

W = round(3 * sigma);
S = exp(-(-W:W).^2 ./ sigma.^2);
S = S / sum(S);

% avoid division by zero in next step
sf(1) = sf(2);
sb(end) = sb(end-1);

% calculate forward and backward deviation
% (convoluted with Gaussian to smoothe signal)
df = (mf - signal)./sf;
db = (signal - mb)./sb; 
d = conv(df + db, S);

% find maximum
[dmax,index] = max(d);

% check against threshold
if nargin < 3
    threshold = 4;
end

if dmax > threshold
    index = index - W;
else
    index = length(signal);
end
