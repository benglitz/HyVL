function OUT = C_parseID(IN)

if ~isempty(IN) && isnumeric(IN{1}) && isnumeric(IN{2})
  ANum = IN{1};
  RNum = IN{2};
  OUT = {'Animal',['mouse',num2str(ANum)],'Recording',RNum}; 
  if length(IN)>2 OUT = [OUT(:)',IN(3:end)]; end
else
  OUT = IN;
end