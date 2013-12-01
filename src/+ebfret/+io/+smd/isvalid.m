function v = isvalid(smd)
% isvalid(smd)
%
% Checks if supplied struct is a valid single-molecule 
% dataset instance
%
% Jan-Willem van de Meent
% $Revision: 1.0$  $Date: 2012/10/16$

v = false;
if ~isstruct(smd)
    return
end
if ~all(isfield(smd, {'type', 'id', 'attr', 'columns', 'data'}))
    return
end
if ~all(isfield(smd.data, {'id', 'attr', 'index', 'values'}))
    return
end
v = true;