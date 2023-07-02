function Vec = C_addRamp(Vec,RampDur,SR,Style,StartLevel)

if nargin < 4; Style = '<=>'; end
if nargin < 5; StartLevel = 0; end

RampSteps=floor(RampDur*SR);
Time = linspace(0,pi,RampSteps);
if size(Vec,1)>size(Vec,2);  Time = Time'; end
  
StartRamp = (1-cos(Time))/2; % Invert cosine to create the ramp at the beginning
if StartLevel; StartRamp = rescale(StartRamp,StartLevel,1); end % rescale ramp if starting signal level specified
EndRamp = StartRamp(end:-1:1); % Invert in time to create the ramp at the end

switch Style
 case '<==' % Apply Ramp only at the beginning
  Vec(1:RampSteps) = StartRamp.*Vec(1:RampSteps);
 case '==>' % Apply Ramp only at the end
  Vec(end-RampSteps+1:end) = EndRamp.*Vec(end-RampSteps+1:end);
 case '<=>' % Apply Ramp at both ends
  Vec(1:RampSteps) = StartRamp.*Vec(1:RampSteps);
  Vec(end-RampSteps+1:end) = EndRamp.*Vec(end-RampSteps+1:end);
  otherwise error('Unknown Rampshape!');
end
