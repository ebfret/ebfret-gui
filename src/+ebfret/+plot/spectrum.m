function map = spectrum(m, varargin)
%map = spectrum(m)
%   Dark rainbow colormap. Returns size [m 3] array of colormap 
%   colors. If m is not specified, the current figure's colormap 
%   length is used. 

if nargin < 1
    m = size(get(gcf,'colormap'),1); 
end

ip = inputParser();
ip.StructExpand = true;
ip.addParamValue('min_hue', 0.3, @isnumeric);       
ip.addParamValue('max_hue', 1.1, @isnumeric);       
ip.addParamValue('min_vel', 0.66, @isnumeric);       
ip.addParamValue('max_vel', 0.66, @isnumeric);       
ip.parse(varargin{:});
args = ip.Results;

h = mod(linspace(args.min_hue, args.max_hue, m), 1);
s = ones(m,1);
v = args.max_vel - (args.max_vel - args.min_vel) .* linspace(-1, 1, m).^2;

if isempty(h)
  map = [];
else
  map = hsv2rgb([h(:) s(:) v(:)]);
end
