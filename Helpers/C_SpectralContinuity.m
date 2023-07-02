function SC = C_SpectralContinuity(S,NoiseThresh)

TShiftMax = 4;
% COMPUTE MAXIMUM ACROSS FREQUENCIES
[MAX,MAXPos] = max(S);
MAX = full(MAX); MAXPos = full(MAXPos);

% COMPUTE DIFFERENCES IN MAXIMUM POSITION BETWEEN TIME STEPS
MAXPosDiff = zeros(size(MAX));
MAXPosDiff(1:end-1) = abs(MAXPos(1:end-1) - MAXPos(2:end));
% IGNORE SMALL MAXIMA = NOISE
MAXPosDiff(MAX(1:end-1)<NoiseThresh) = size(S,1);

% INTEGRATE DIFFERENCES ALONG TIME IN BOTH DIRECTIONS
CSFwd = cumsum(MAXPosDiff); 
CSBwd = cumsum(MAXPosDiff(end:-1:1)); 

% COMPUTE TOTAL DIFFERENCE AT TShiftMax in both directions
SCFwd = CSFwd(TShiftMax+1:end)-CSFwd(1:end-TShiftMax);
SCBwd = CSBwd(TShiftMax+1:end)-CSBwd(1:end-TShiftMax);
SC = size(S,1)*ones(size(MAX));

% ASSIGN THE FORWARD INTEGRATION 
SC(1:end-TShiftMax) = SCFwd;
% COMPARE WITH THE BACKWARD INTEGRATION, KEEP LOWER VALUE
SC(TShiftMax+1:end) = min(SC(TShiftMax+1:end),SCBwd(end:-1:1));
