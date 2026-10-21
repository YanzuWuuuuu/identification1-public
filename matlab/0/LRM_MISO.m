function varargout=LRM_MISO(varargin)
% Compatibility entry point for the shared, scaled local rational fit.
previousPath=path; cleanup=onCleanup(@()path(previousPath));
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'common'));
[varargout{1:nargout}]=lrm_miso_core(varargin{:});
end
