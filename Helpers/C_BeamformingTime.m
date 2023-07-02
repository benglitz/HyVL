function ResAll = C_Beamforming(varargin)

P = parsePairs(varargin);
checkField(P,'Data'); % as Time X Microphones
checkField(P,'MicrophonePositions'); % as NMic X 3, in absolute, setup coordinates
checkField(P,'SR'); % Sampling Rate
% Scanning parameters
checkField(P,'Method','DelaySum');
checkField(P,'FRange',[50000,80000]); % Scanning frequency for the beamforming
checkField(P,'XRange',[-0.2,0.2]); % In absolute, setup coordinates
checkField(P,'YRange',[-0.15,0.15]); % In absolute, setup coordinates
checkField(P,'ZDist');
checkField(P,'XStep',0.01);
checkField(P,'YStep',0.01);

% Analysis Parameters
checkField(P,'TBin',0.02); % Analysis Bin in Seconds
checkField(P,'TStep',P.TBin);
checkField(P,'SpeedOfSound',343); % m/s
checkField(P,'FIG',1);
checkField(P,'TRange',[])
checkField(P,'Verbose',0)
checkField(P)

[NSteps,NMic] = size(P.Data);

% ASSIGN MICROPHONE POSITIONS
Xm = P.MicrophonePositions(:,1);
Ym = P.MicrophonePositions(:,2);
Zm = P.MicrophonePositions(:,3);

% ASSIGN BEAMFORMING LOCATIONS
X = P.XRange(1):P.XStep:P.XRange(2);  NStepsX = length(X);
Y = P.YRange(1):P.YStep:P.YRange(2);  NStepsY = length(Y);
ResAll.x = X; ResAll.y = Y;

% LOOP OVER TIME TO CREATE BEAMFORMING ESTIMATE
Duration = NSteps/P.SR;
if ~isempty(P.TRange)
  if length(P.TRange) == 1;     TimeEst = P.TRange;
  else                                          TimeEst = [P.TRange(1):P.TStep:P.TRange(2)];
  end
else TimeEst = 0:P.TStep:Duration-P.TBin;
end
NTimeEst = length(TimeEst);
NStepsBin = round(P.TBin*P.SR);
ResAll.Time = TimeEst;

% BANDPASS FILTER THE DATA
if ~isempty(P.FRange)
  [bBP,aBP] = fir1(128,P.FRange([1,end])/P.SR*2);
  P.Data = filter(bBP,aBP,P.Data);
end

if P.Verbose; fprintf(['Starting estimate for ',num2str(NTimeEst),'\n']); end
for iT = 1:(NTimeEst)
  if mod(iT,10)==0; fprintf([num2str(iT),' ']); end
  iStart = round(TimeEst(iT)*P.SR)+1;
  iStop = iStart + NStepsBin-1;
  cInd  = [iStart:iStop];
  cData = double(P.Data(cInd,:));

  switch P.Method
    case 'DelaySum'
      % Compute Delay and Sum nearfield solution
      Res = nearfieldDelayAndSum(cData,Xm,Ym,Zm,X,Y,P.ZDist,P.SpeedOfSound,P.SR);
    case 'DelaySumGPU'
      % Compute Delay and Sum nearfield solution on the GPU
      Res = nearfieldDelayAndSumGPU2(cData,Xm,Ym,Zm,X,Y,P.ZDist,P.SpeedOfSound,P.SR);
    otherwise
      error('Beamforming technique not implemented.');
  end
  ResAll.ProjectionXY(:,:,iT)  = Res'; % Transposed since Cam64 Portal Data also in YX and does not require transpose for plotting

end
