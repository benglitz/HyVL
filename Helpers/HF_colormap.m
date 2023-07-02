function M = HF_colormap(Colors,Division,N)

if ~exist('N','var') | isempty(N) N = 255; end
if ~exist('Division','var') | isempty(Division) Division = [0:length(Colors)-1]; end

Ncolors = length(Colors); ColorSteps = diff(Division); 
ColorSteps = [0,round(N*cumsum(ColorSteps)/sum(ColorSteps))];
M = zeros(N,3);

for i=2:Ncolors
  ColStart = ColorSteps(i-1)+1;
  ColStop = ColorSteps(i);
  for j=1:3
    M(ColStart:ColStop,j) ...
      = linspace(Colors{i-1}(j),Colors{i}(j),ColorSteps(i)-ColorSteps(i-1));
  end  
end