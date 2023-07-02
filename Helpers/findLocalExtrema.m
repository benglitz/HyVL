function [pos,vals] = findLocalExtrema(Vec,Opt,N)
% findLocalExtrema: looks for local extrema 
% Opt : can be either 'min' or 'max', delivering the respective type of extremum
% Boundary Extrema are left out, since we do not have information on the 
% behavior past the boundarz
pos = []; k=2;
dsdVec = diff(sign(diff(Vec)));
if Opt == 'max'
  pos = find(dsdVec==-2)+1;
elseif Opt=='min'
  pos = find(dsdVec== 2)+1;
else error('unknown option Opt!');
end
pos=pos(pos<length(Vec));
vals = Vec(pos);
if exist('N','var') & ~isempty(pos) pos = pos(1:N); vals = vals(1:N); end
