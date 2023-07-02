function printupdate(varargin)

switch nargin
  case 3; Mode = 'Numbers';
    S = sprintf('%d / %d',varargin{1},varargin{2}); Init = varargin{3};
  case 2; Mode = 'String';
    S = varargin{1}; Init = varargin{2};
  case 1; Mode = 'String'; S = varargin{1}; Init = 0;
end

persistent NWritten
if Init NWritten = 0; end
if ~isempty(NWritten)
  fprintf(repmat('\b',1,NWritten));
end
NWritten = fprintf(S);
