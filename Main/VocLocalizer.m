function R = VocLocalizer(varargin)
% Localize vocalizations in space using a microphone array
%
% Notes for 1D Measurements:
% Realignment: 29.9.2016 (Jesse/Bernhard)
% - We think the intended distance between the lower edges of the microphones was
%   - Horizontal: 46cm, i.e. 23cm from the center (this corresponds to the inner edges of the rear platform walls) 
%   - Vertical :     15cm above the upper end of the platform 
%                      = 34.4cm above the platform level
%   - Depth : centered over the platforms
%   - Microphone angle : 45 degrees
%   NOTE: This means that the membranes of the microphones are at a
%   different location. Since the diameter sof the microphone front is 36mm,
%   the center of the membrane is located 1.27cm further up and inward.
%  This means the position of the MEMBRANES is 
%  - Horizontal: 43.46cm, i.e. 21.73cm from the center (this corresponds to the inner edges of the rear platform walls) 
%   - Vertical :     16.27cm above the upper end of the platform 
%                      = 35.67cm above the platform level
%   - Depth : centered over the platforms
%   - Membrane angle : 45 degreess
%  The Membrane position should be the authorative position of the
%  microphones.
% 
% The Height of the speaker in relation to the microphones was remeasured as 31.4cm
%
% Old Notes on Position:
% MicDist = 0.46; %m % Measurements from new setup
% - Estimate of previous distances : 0.38
% MicHeight = 0.314; %m % Measurement of Speaker
% MicHeight = 0.334; %m % Measurement for Mouse
% - Estimate of previous height: 0.35
% Platform is assumed to be at Z = 0.
%
% Lenses Used:
% - 12.5mm : leads to 16.4 mm % measured on 20.1.2017
% - 35 mm (XR300, SainSonic) : leads to 56mm 
%
% TODO: 
% - Test combination of Specgram and EWGCC for cases, whenever the localization is close to the center between the microphones

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'Sounds')
checkField(P,'SR')
checkField(P,'NFFT',250);
checkField(P,'RampDur',0.002); % 2ms ramp
checkField(P,'HighPass',35000);
checkField(P,'LowPass',120000);
checkField(P,'RemoveAmplitude',0);
checkField(P,'RemoveMean',1); % New Option!
checkField(P,'CorrMethod','GCC');
checkField(P,'Offset',1);
checkField(P,'MaxMethod','Max');
checkField(P,'CorrTime',[]); % 0.1 ms
checkField(P,'DelayRange',P.CorrTime); %75us
checkField(P,'MicrophonePositions',{[-0.23,0,0.314],[0.23,0,0.314]})
checkField(P,'SourceHeight',0);
checkField(P,'EstimationMethod',[]); % Analytical or General
checkField(P,'Estimation',[]); % X , Y , or XY
checkField(P,'AnalyticalMethod','Geometric');
checkField(P,'CenterShift',0);
checkField(P,'GainCorrection',1.2);
checkField(P,'IntersectionMethod','Lines');
checkField(P,'VSound',343); %m/s
checkField(P,'StimulusPosition',[]); % For checking in simulated recordings
checkField(P,'DelayJitter',0); % Introduce artificial jitter, only used for illustration in methods figure
checkField(P,'Filter',[])
checkField(P,'FIG',1); %
checkField(P);

%% PREPARE ANALYSIS
NMic = length(P.MicrophonePositions);
NMicPairs = NMic*(NMic-1)/2;
for iM=1:NMic; 
  P.Sounds{iM} = vertical(P.Sounds{iM}); 
  P.MicrophonePositions{iM}(3) = P.MicrophonePositions{iM}(3) - P.SourceHeight; 
end
NSteps = length(P.Sounds{1});
if ischar(P.CorrMethod) P.CorrMethod = {P.CorrMethod}; end

if isempty(P.CorrTime)
  k=0; InterMicDist = zeros(NMicPairs,1);
  for iM1=1:NMic-1; for iM2=iM1+1:NMic 
      k=k+1;
      InterMicDist(k) = norm(P.MicrophonePositions{iM1} - P.MicrophonePositions{iM2}); 
    end; end
  P.CorrTime = 1.3*max(InterMicDist(:)/2/P.VSound); 
  % Division by two since it is measured from the middle, multiplication by 1.3, since the delay can be bigger due to the Z-difference
end
if isempty(P.DelayRange);  P.DelayRange = P.CorrTime;
else  P.DelayRange = min(P.DelayRange,P.CorrTime); 
end

CorrSteps = round(P.CorrTime*P.SR);

if P.RemoveMean
  for iS = 1:length(P.Sounds)
    P.Sounds{iS} = P.Sounds{iS} - mean(P.Sounds{iS});
  end
end

%% ADD GATING RAMPS
if P.RampDur
  for iS=1:length(P.Sounds)
    P.Sounds{iS} = C_addRamp(P.Sounds{iS},0.002,P.SR,'<=>');
  end
end

%% POTENTIALLY FILTER THE DATA
%warning('Introduce a bandpass filtfilt here before continuing based on FMin/FMax');
if P.HighPass
  [b,a]=butter(6,P.HighPass/(P.SR/2),'high');
  for i=1:length(P.Sounds);  P.Sounds{i} = filter(b,a,P.Sounds{i}); end
end
if P.LowPass
  [b,a]=butter(6,P.LowPass/(P.SR/2),'low');
  for i=1:length(P.Sounds);  P.Sounds{i} = filter(b,a,P.Sounds{i}); end
end
% if isempty(P.Filter) % Removed: worked well if the filter is 
%   warning('Filter not set externally');
%   P.Filter = designfilt('bandpassiir','FilterOrder',20, 'HalfPowerFrequency1',30000,'HalfPowerFrequency2',100000, 'SampleRate',P.SR);
% end
% for i=1:length(P.Sounds);  P.Sounds{i} = filter(P.Filter,P.Sounds{i}); end

if P.RemoveAmplitude
  for i=1:length(P.Sounds)
    H{i} = abs(hilbert(P.Sounds{i})); P.Sounds{i} = P.Sounds{i}./H{i};
  end
end

% PREPARE VARIABLES
SignalToNoiseAll = NaN*ones(NMic-1,NMic);

%% COMPUTE CENTER DISTANCE FOR EACH PAIR OF MICROPHONES
P.PlotToFigure = 0;
if isnumeric(P.FIG) & mod(P.FIG,1)==0 & P.FIG>0 P.PlotToFigure =1; figure(P.FIG); clf; set(P.FIG,'Position',[300,40,1000,400]); end 
if P.PlotToFigure
  DC = axesDivide(2,1);
  [~,P.AH] = axesDivide(1,NMicPairs,DC{1},[],[0.3],'c');
  linkaxes(P.AH,'y');
  P.AH(end+1) = axes('Pos',DC{2});
end
  
MicPairs = {}; iL=0; MAX = 0; % LOOP OVER MICROPHONE PAIRS
for iM1 = 1:NMic
  for iM2 = iM1+1:NMic
    iL = iL+1;
    MicPairs{iL} = [iM1,iM2];
    Sound1 = single(P.Sounds{iM1});
    Sound2 = single(P.Sounds{iM2});
    if P.DelayJitter
      rand('seed',18);
      cJitter = round(P.DelayJitter*rand(2,1))+1;
      Sound1 = Sound1(cJitter(1):end);
      Sound2 = Sound2(cJitter(2):end);
    end
    SCorrFinal = ones(2*CorrSteps+1,1);
    for iM = 1:length(P.CorrMethod)
      switch P.CorrMethod{iM}
        case 'Specgram';
          % COMPUTE SPECTROGRAM OF THE SOUND FROM TWO MICROPHONES
          S1 = HF_specgram(Sound1,P.NFFT,P.SR,[],P.NFFT-P.Offset,0,0);
          S2 = HF_specgram(Sound2,P.NFFT,P.SR,[],P.NFFT-P.Offset,0,0);
          
          % SELECT FREQUENCY INDICES TO COMPUTE THE CROSSCORRELATION
          FreqMarginals = mean(abs(S1)+abs(S2),2);
          FreqInds = find( FreqMarginals > (min(FreqMarginals)+max(FreqMarginals)/2));
          
          % COMPUTE CROSSCORRELATIONS
          Corrs= zeros(length(FreqInds),2*CorrSteps+1);
          for iF=1:length(FreqInds)
            cInd = FreqInds(iF);
            Corrs(iF,:) = xcorr(abs(S1(cInd,:)),abs(S2(cInd,:)),CorrSteps,'unbiased');
          end
          Corrs = Corrs./repmat(FreqMarginals(FreqInds),1,size(Corrs,2));
          
          XCorr  = mean(Corrs)';
          DT = [-CorrSteps:CorrSteps]/P.SR;
          %CorrKernel = 0.01*exp(-DT.^2/(2*0.00065^2)) + 1;
          XCorr = XCorr / max(XCorr);
          %XCorr = XCorr.*CorrKernel';          
          SCorr{iM} =XCorr.^10;
          
        case 'GCC';
          NFFT = NSteps -1 + mod(NSteps,2);
          S1F = fft(double(Sound1),NFFT); S2F = fft(double(Sound2),NFFT);
          F = [1:NFFT]/NFFT*P.SR; F = F(1:round(NFFT/2));
          Equalizer = S1F.*conj(S2F)./(abs(S1F).*abs(conj(S2F)));
          XCorr = abs( ifft(Equalizer) );
          XCorr = [ XCorr(round(NFFT/2)+1:end) ; XCorr(1:round(NFFT/2))  ];
          %XCorr = abs( hilbert( XCorr ) );
          XCorr = XCorr(round(end/2)-CorrSteps:round(end/2)+CorrSteps);
          SCorr{iM} = XCorr;
                    
        case 'EWGCC'
          NFFT = NSteps -1 + mod(NSteps,2);
          S1F = fft(double(Sound1),NFFT); S2F = fft(double(Sound2),NFFT);
          F = [1:NFFT]/NFFT*P.SR; F = F(1:round(NFFT/2));
          Step = 10; FR = F(Step/2:Step:end);
          AThresh = 0.018;
          NBins = floor(NFFT/Step);
          NStepsBinned = NBins*Step;
          A1R = reshape(abs(S1F(1:NStepsBinned)),Step,NBins); mA1R = mean(A1R);
          A2R = reshape(abs(S2F(1:NStepsBinned)),Step,NBins); mA2R = mean(A2R);
          mA1 = interp1(FR,mA1R(1:length(FR)),F);
          mA2 = interp1(FR,mA2R(1:length(FR)),F);
          SmallInd = (mA1<AThresh*max(mA1)) | (mA2<AThresh*max(mA2));
          SmallInd = [SmallInd,flipud(SmallInd)];
          Equalizer = S1F.*conj(S2F)./(abs(S1F).*abs(conj(S2F)));
          Equalizer(SmallInd) = 0;
          XCorr = abs(ifft(Equalizer));
          XCorr = [ XCorr(round(NFFT/2)+1:end) ; XCorr(1:round(NFFT/2))  ];
          XCorr = XCorr(round(end/2)-CorrSteps:round(end/2)+CorrSteps);
          %HXCorr = abs( hilbert( XCorr ) );
          XCorr = LF_removeGaussian(XCorr,1.5e-4,P.SR);
          SCorr{iM} = XCorr;
          
        case 'XCorr'
          SCorr{iM} = xcorr(Sound1,Sound2,CorrSteps,'unbiased');
          SCorr{iM} = abs(hilbert(SCorr{iM}));
          
        case 'Micro'
          OversampleFactor = 8;
          % OVERSAMPLE DATA
          SoundsLong{1} = interpft(Sound1,NSteps*OversampleFactor+1);
          SoundsLong{2} = interpft(Sound2,NSteps*OversampleFactor+1);
          P.SR = P.SR*OversampleFactor;
          
          % REMOVE PHASE
          for i=1:length(SoundsLong)
            H{i} = abs(hilbert(SoundsLong{i}));
            SoundsLong{i} = SoundsLong{i}./H{i};
          end
          
          % CONSTRUCT HISTOGRAM OF THE PHASE MATCHES
          CorrSteps = round(P.CorrTime*P.SR/2);
          WindowSteps = P.CorrTime*P.SR;
          NFFT = WindowSteps;
          MaxHist = zeros(2*CorrSteps+1,1);
          Steps = round(linspace(1,NSteps-WindowSteps,100));
          for i=1:length(Steps)
            cStep = Steps(i);
            cInd = cStep+[0:WindowSteps-1];
            S1 = SoundsLong{1}(cInd); S2 = SoundsLong{2}(cInd);
            X = xcorr(S1,S2,CorrSteps,'biased');
            X(X<max(X)/2) = -inf;
            [StimPosMean,Vals] = findLocalExtrema(X,'max');
            MaxHist(StimPosMean) = MaxHist(StimPosMean) + 1;
          end
          SCorr{iM} = MaxHist;
          
        otherwise error(['Method ',P.CorrMethod,' not known.']);
      end
      
      % COMBINE THE ESTIMATES FROM DIFFERENT METHODS
      SCorrFinal = SCorr{iM}.*SCorrFinal;
    end
    
    %% POST PROCESS LOCALIZATION
    % SMOOTH THE CORRELATION
    %SCorrFinal = relaxc(SCorrFinal,2);
   % SCorrFinal = SCorrFinal-min(SCorrFinal); SCorrFinal = SCorrFinal/max(SCorrFinal);
    SCorrAll{iM1,iM2} = SCorrFinal;
    
    % SUBSELECT A RANGE FOR CHOOSING THE RIGHT DELAY
    CorrTime = [-CorrSteps:CorrSteps]/P.SR;
    NotInd = find(abs(CorrTime)>P.DelayRange);
    SCorrRange = SCorrFinal;
    SCorrRange(NotInd) = -inf;
    
    if P.PlotToFigure
      axes(P.AH(iL)); MAX = max([MAX, max(SCorrFinal)]);
      plot(CorrTime,SCorrFinal,'r');  %plot(CorrTime,XCorr,'k');
      title(['Mic. ',num2str(iM1),' vs. Mic. ',num2str(iM2)]);
      ylim([0,MAX]);
    end
    
    % EXTRACT CENTER OF CORRELATION
    switch P.MaxMethod
      case 'Weighted';
        M = max(SCorrRange);
        SCorrRange(SCorrRange<0.75*M) = 0;
        StimPosMean = sum([1:length(SCorrRange)]'.*SCorrRange./(sum(SCorrRange)));
      case 'Max';
        [M,StimPosMean] = max(SCorrRange);
      otherwise error('Maximum Method not known');
    end
    SignalToNoiseAll(iM1,iM2) = (M - std(SCorrRange(~isinf(SCorrRange)))) / std(SCorrRange(~isinf(SCorrRange)));
    
    % CONVERT TIME TO POSITION ON CONNECTING LINE
    DeltaTime(iM1,iM2) = real((StimPosMean-(CorrSteps+1))/P.SR); % DeltaTime is negative, if signal arrives first on the left (channel 1), i.e. position is also on the left
    DeltaDist(iM1,iM2) = P.VSound * DeltaTime(iM1,iM2);
    MidDist(iM1,iM2) = DeltaDist(iM1,iM2)/2;

    if abs(imag(MidDist(iM1,iM2)))>0
      MidDist(iM1,iM2) = NaN; fprintf('Warning: Out of bounds distance!\n');
    end
        
  end
end


R.SignalToNoiseAll = SignalToNoiseAll;
R.MidDist = MidDist;
R.SCorrAll = SCorrAll;
R.CorrTime = CorrTime;
R.DeltaTime = DeltaTime;


%% TRANSLATE POSITIONS ALONG THE MICROPHONE-CONNECTING LINES
% TO OFFCENTER POSITIONS
for iM=1:length(P.MicrophonePositions)
  MicPos(iM,:) = P.MicrophonePositions{iM};
end
R.MIN = min(MicPos) - 0.05; R.MAX = max(MicPos) + 0.05;

if isempty(P.Estimation)
  if length(unique(MicPos(:,1))==1); P.Estimation = 'Y'; end;
  if length(unique(MicPos(:,2))==1); P.Estimation = 'X'; end;
  if length(unique(MicPos(:,1)))>1 && length(unique(MicPos(:,2)))>1; P.Estimation = 'XY'; end
end
if isempty(P.EstimationMethod)
  if length(P.Estimation) == 1 % only X or Y Estimation
    P.EstimationMethod = 'Analytical';
  else  
    P.EstimationMethod = 'General';
  end
end

switch P.EstimationMethod
  case 'Analytical';
    R = LF_analyticalEstimation(DeltaDist(1,2),R,P);
    
  case 'General';
    R = LF_generalEstimation(MidDist,MicPos,MicPairs,R,P);
end

%% COMPARE TO REAL POSITION
if ~isempty(P.StimulusPosition)
  Error = norm(R.StimPosMean - P.StimulusPosition(1:2));
  fprintf(['Error = ',num2str(Error*1000),'mm (Estimated: ',num2str(R.StimPosMean),'mm)\n']);
end

%% PLOT RESULTS
if ~isempty(P.FIG) && P.FIG ~= 0 
  R.Sounds = P.Sounds;
  R.Time = [1:NSteps]/P.SR;
  R.SCorr = SCorrAll;
  LF_showLocalization(P.StimulusPosition,R.StimPosMean,P.MicrophonePositions,R,P); 
end

%CorrMidDist{iM1,iM2} = LF_translateTime2Space(CorrTime,P.LocMethod,P.CenterShift,P.GainCorrection,P.MicrophonePositions([iM1,iM2]));
%R.CorrMidDist = CorrMidDist;

end % END OF MAIN FUNCTION

%% COMPUTE POSITION FROM TIMING
% ANALYTICAL ESTIMATION
function R = LF_analyticalEstimation(DeltaDist,R,P)

switch P.AnalyticalMethod
  case 'Geometric';
    % COMPUTE POSITION BASED ON SETUP GEOMETRY AND CORRELATION
    if P.MicrophonePositions{1}(3) ~= P.MicrophonePositions{2}(3); error('Microphones need to be at the same height'); end 
    MicHeight = P.MicrophonePositions{1}(3);
    MicDist = P.MicrophonePositions{2}(1)-P.MicrophonePositions{1}(1);
    T = DeltaDist;   H = MicHeight;  D  = MicDist;
    
    StimPosMean = 0.5.*T.* sqrt((4*H^2 + D^2 - T.^2)./(D^2-T.^2));   
    StimPosMean = StimPosMean - P.CenterShift;
     
  case 'Empirical';
    StimPosMean = (DeltaDist - P.CenterShift)/P.GainCorrection;

  otherwise error(['Localization Method ',LocMethod,'not known.']);
end
if imag(StimPosMean(1)) keyboard; end
R.StimPosMean = [StimPosMean,P.MicrophonePositions{1}(2)];
R.StimPosSD = [NaN,NaN];
end

%% GENERAL ESTIMATION
function R = LF_generalEstimation(MidDist,MicPos,MicPairs,R,P)
% MicSep : Separation between the microphones (D in the manuscript)
% OffCenter : Vertical Distance between the sound source and microphone membrane (H in the manuscript)
%                    below inserted as Platform dist.
% MidDist : Difference in PathLength, converted via the speed of sound (Delta P in the manuscript)
InvFun = @(MidDist,OffCenter,MicSep)MidDist * sqrt((4*(OffCenter).^2 + MicSep.^2 - MidDist.^2)./(MicSep.^2 - MidDist.^2));
Range = sqrt(2)*max(R.MAX-R.MIN);
OrthDist = [-Range:0.001:Range];
MicHeight = unique(MicPos(:,3)); assert(length(MicHeight)==1,'Microphones are assumed to be at the same height');
PlatformDist = sqrt(OrthDist.^2 + MicHeight^2); % Absolute, orthogonal distance to the shifted line of potential positions

%% COMPUTE THE SPATIAL LOCATION FROM THE VARIOUS MEASUREMENTS
% Compute the lines defined by each Microphone Pair
for iP = 1:length(MicPairs)
  P1 = P.MicrophonePositions{MicPairs{iP}(1)}(1:2);
  P2 = P.MicrophonePositions{MicPairs{iP}(2)}(1:2);
  % Find Orthogonal Vector onto Line between Microphones
  VA = (P2-P1); VA = VA/norm(VA);
  VAO = [VA(2),-VA(1)]; VAO = VAO/norm(VAO);
  % Set Starting Point based on MidDist
  cMidDist = MidDist(MicPairs{iP}(1),MicPairs{iP}(2));
  Lines(iP).P0 = (P1 + P2)/2;
  Lines(iP).P1 = P1;
  Lines(iP).P2 = P2;
  Lines(iP).V = VAO;
  Lines(iP).Mics = MicPairs{iP};
  Lines(iP).MidDist = cMidDist;
  Lines(iP).MicSep = norm(P2-P1);
    
  % COMPUTE THE MID
  Lines(iP).MidDists = real(InvFun(Lines(iP).MidDist,PlatformDist,Lines(iP).MicSep));
  
  % ROTATION MATRIX for each pair of microphones
  % using the orthogonal vector onto the connecting line
  alpha = angle(Lines(iP).V(1) + sqrt(-1)*Lines(iP).V(2));
  Rotator = real([cos(alpha),-sin(alpha);sin(alpha),cos(alpha)]);
  Lines(iP).MidDistsRot = Rotator*[OrthDist;Lines(iP).MidDists] + repmat([Lines(iP).P0(1);Lines(iP).P0(2)],[1,length(OrthDist)]);
end
R.Lines = Lines;

%% ESTIMATE OVERLAP OF ESTIMATION LINES
switch P.Estimation
  case 'XY'; NStepsX = 1001; NStepsY=1001;
  case 'X'; NStepsX = 5001; NStepsY = 11; % MIN(2)=-1e-3; MAX(2) = 1e-3;
  case 'Y'; NStepsY = 5001; NStepsX = 11; % MIN(1)=-1e-3; MAX(1) = 1e-3;
end
mX = linspace(R.MIN(1),R.MAX(1),NStepsX);
mY = linspace(R.MIN(2),R.MAX(2),NStepsY);
M = zeros(length(mY),length(mX));
% FILL A MATRIX WITH THE LINES LOCALLY PLACED
for iP=1:length(MicPairs)
  cP = Lines(iP).MidDistsRot;
  i1 = ceil((cP(2,:)-R.MIN(2))/(R.MAX(2)-R.MIN(2))*NStepsY);
  i2 = ceil((cP(1,:)-R.MIN(1))/(R.MAX(1)-R.MIN(1))*NStepsX);
  SelInd = i1>0 & i1<=NStepsY & i2>0 & i2<=NStepsX;
  i1 = i1(SelInd); i2 = i2(SelInd);
  Ind = sub2ind(size(M),i1,i2);
  M(Ind)  = M(Ind) + 1;
end
Steps = [-2:0.1:2];
[X,Y] = meshgrid(Steps,Steps);
K = exp(-sqrt(X.^2+Y.^2)/0.6);
K = K/(1.1*norm(K(:),2));
M = filter2(K,M)/length(MicPairs); % Normalize by the maximally possible height of the Interasection
R.M = M; R.mX = mX; R.mY = mY;

% IntersectionMethod for: 
% - Platform2D (3 Microphones): Intersections
% - Platform2DFoxp2 (4 Microphones): Map

IntersectionMethod = 'Map';
if strcmp(IntersectionMethod,'Intersections')
  % Compute the intersection points between all pairs of lines
  for iP1=1:length(MicPairs)-1
    for iP2=iP1+1:length(MicPairs)
      C1 = Lines(iP1).MidDistsRot;
      C2 = Lines(iP2).MidDistsRot;
      [cInterX,cInterY] = intersections(C1(1,:),C1(2,:),C2(1,:),C2(2,:));
      if length(cInterX)>1
        for iI = 1:length(cInterX) Norms(iI) = sqrt(cInterX(iI).^2+cInterY(iI).^2); end
        [~,cInd] = min(Norms);
        cInterX = cInterX(cInd); cInterY = cInterY(cInd);
      end
      R.Intersections{iP1,iP2} = [cInterX,cInterY]';
      R.IntersectionsByMic{iP1,iP2}{1} = Lines(iP1).Mics;
      R.IntersectionsByMic{iP1,iP2}{2} = Lines(iP2).Mics;
    end
  end
end

%% EXTRACT PROJECTIONS ONTO THE ESTIMATED DIRECTIONS
switch P.Estimation
  case 'X'
    iY = find(mY==min(abs(mY)));
    iX = find(M(iY,:)==max(M(iY,:)));
    R.StimPosMean = [mean(mX(iX)),MicPos(1,2)]; R.StimPosSD = [std(mY(iY)),0];
  case 'Y'
    iX = find(mX==min(abs(mX)));
    iY = find(M(iX,:)==max(M(iX,:)));
    R.StimPosMean = [MicPos(1,1),mean(mY(iY))];  R.StimPosSD = [0,std(mY(iY))];
  case 'XY'
    switch IntersectionMethod
      case 'Map'
        [iY,iX] = find(M>0.9*max(M(:)));
        R.StimPosMean = [mean(mX(iX)),mean(mY(iY))];
        R.StimPosSD     = [std(mX(iX)),std(mY(iY))];
      case 'Intersections'
        % Find the two microphone pairs with the highest signal to noise in their correlation
        Values = sort(R.SignalToNoiseAll(:),'descend'); Values = Values(~isnan(Values));
        [M1(1),M1(2)] = find(R.SignalToNoiseAll == Values(1));
        [M2(1),M2(2)] = find(R.SignalToNoiseAll == Values(2));
        for iP1=1:size(R.Intersections,1)
          for iP2=1:size(R.Intersections,2)
            if ~isempty(R.IntersectionsByMic{iP1,iP2}) & ...
                ((R.IntersectionsByMic{iP1,iP2}{1} == M1 & R.IntersectionsByMic{iP1,iP2}{2} == M2) ...
                | (R.IntersectionsByMic{iP1,iP2}{1} == M2 & R.IntersectionsByMic{iP1,iP2}{2} == M1))
              R.StimPosMean = R.Intersections{iP1,iP2}';
            end
          end
        end
        if isempty(R.StimPosMean);  R.StimPosMean = mean([R.Intersections{:}],2)'; end
        R.StimPosSD = std([R.Intersections{:}],[],2)';
    end
end
end

%% SHOW THE LOCALIZATION RESULTS
function Pos = LF_showLocalization(StimPos,StimPosMean,MicPos,R,P)

if P.PlotToFigure; AH = P.AH(end); else AH = P.FIG; end
axes(AH); hold on;

% PLOT MICROPHONES & STIMULUS
NMic = length(MicPos);
for iM = 1:NMic
  plot3(MicPos{iM}(1),MicPos{iM}(2),MicPos{iM}(3),'.k','MarkerSize',15,'HIttest','off');
  text(MicPos{iM}(1),MicPos{iM}(2),1.5*MicPos{iM}(3),['M',num2str(iM)],...
    'ButtonDownFcn',{@LF_showData,R.Time,R.Sounds{iM}},'Horiz','center','FontSize',7);
  plot3([MicPos{iM}(1),MicPos{iM}(1)],[MicPos{iM}(2),MicPos{iM}(2)],[0,MicPos{iM}(3)],'--','Color',[0.5,0.5,0.5],'HIttest','off');
end
if ~isempty(StimPos)
  plot3(StimPos(1),StimPos(2),StimPos(3),'o','MarkerSize',6,'Hittest','off','Color',[0.5,0.5,1]);
end
axis square; hold on; box on; grid on;
axis([R.MIN(1),R.MAX(1),R.MIN(2),R.MAX(2),0,R.MAX(3)]) ;

% PLOT THE INTERSECTION LINES
if isfield(R,'Lines')
  for iP = 1:length(R.Lines)
    cLine = R.Lines(iP);
    %plot3([cLine.P1(1),cLine.P2(1)],[cLine.P1(2),cLine.P2(2)],[0,0],'--','Color',[0.5,0.5,0.5],'HIttest','off');
    plot3([cLine.P1(1),cLine.P2(1)],[cLine.P1(2),cLine.P2(2)],...
      [MicPos{cLine.Mics(1)}(3),MicPos{cLine.Mics(2)}(3)],'--','Color',[0.5,0.5,0.5],'HIttest','off');
    %cP = cLine.MidDistsRot;
    %plot(cP(1,:),cP(2,:),'Color',[1,0,0],'Hittest','off');
   % text(cP(1,round(end/2)),cP(2,round(end/2)),['X_{',num2str(cLine.Mics(1)),num2str(cLine.Mics(2)),'}'],...
   %   'ButtonDownFcn',{@LF_showData,R.CorrTime,R.SCorr{cLine.Mics(1),cLine.Mics(2)}},'FontSize',7);
  end
  
  % PLOT THE INTERSECTION POINTS
  if isfield(R,'Intersections')
    for iP1=1:length(R.Lines)-1
      for iP2=iP1+1:length(R.Lines)
        if ~isempty(R.Intersections{iP1,iP2})
          plot3(R.Intersections{iP1,iP2}(1),R.Intersections{iP1,iP2}(2),0,'g.','MarkerSize',10);
        end
      end
    end
  end    
  % PLOT THE POSITION ESTIMATION BASIS
  switch P.IntersectionMethod
    case 'Lines';
      imagesc(R.mX,R.mY,R.M,'Hittest','off');
    case 'Intersections';
  end
end

% PLOT ESTIMATED STIMULUS
plot3(StimPosMean(1),StimPosMean(2),0.000001,'.k','MarkerSize',4);

colormap(gca,HF_colormap({[1,1,1],[1,0,0]}));
switch P.Estimation;
  case 'XY'; view(20,20);
  otherwise view(0,90);
end
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
end

function X = LF_removeGaussian(X,Width,SR)

CorrSteps = floor(length(X)/2);
T = [-CorrSteps:CorrSteps]'/SR;
G = exp(-T.^2/(2*Width.^2));
MeanRange = round(Width*SR);
Ind = round(MeanRange/2):MeanRange;
A = mean(X([-Ind,Ind]+CorrSteps+1));
B = mean(X([1:50,end-50:end]));
X = X - 1.2*(A-B)*G;
end

%% CALLBACK TO PLOT THE SOUNDS / CROSSCORRELATIONS
function LF_showData(H,E,X,Y)

SelType = get(gcf,'SelectionType');
switch SelType
  case 'normal'; FIG = 1000; figure(FIG);
  case {'alt','extend'}; FIG = round(1e6*rand); figure(FIG)
end

% SHOW RASTER PLOT
plot(X,Y);
end
