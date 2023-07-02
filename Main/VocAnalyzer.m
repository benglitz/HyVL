function Vocs = VocAnalyzer(Vocs,varargin)
%% TODO
% COLLECT PROPERTIES AND PERFORM INITIAL ANALYSIS OF VOCALIZATIONS  
% TODO: 
% - USE SPECTRAL PURITY TO DEFINE START AND STOP OF VOCALIZATION
% - ASSIGN NAN TO THE SPECTRAL LINE FOR ALL POINTS WITH VALUE BELOW A
% CERTAIN LEVEL... these can later be replaced by 0, if necessary
% Alternative strategy: maybe start from the center (highest spectral
% purity/local amplitude, and then work towards outside? 
% E.g. include a set of anchor points, and then 

%% PARSE PARAMETERS
P = parsePairs(varargin);
checkField(P,'WithinPhrasePause',0.3)
checkField(P,'SpectrogramFilter','CurvatureFlow')
checkField(P,'FreqRange',[30000,inf])
checkField(P,'Selection',[1:length(Vocs)])
checkField(P,'PitchSalience',0)
checkField(P,'BorderBins',4)
checkField(P,'EdgeThreshold',0.2)
checkField(P,'ExcludeNonUSVs',0)
% COMPUTING LOCATION
checkField(P,'ComputeLocation',1)
checkField(P,'EstimationMethod','General')
checkField(P,'MicrophonePositions',[])
checkField(P,'CorrRange',[])
checkField(P,'SourceHeight',0)
checkField(P,'CenterShift',0)
%  OUTPUT
checkField(P,'Selector',{})
% CHECKING
checkField(P,'CheckProperties',0);
checkField(P,'Verbose',1);
checkField(P);

%% PREPARE VARIABLES BEFOE THE LOOP
Vocs = vertical(Vocs);
IVIs = [Vocs.Interval];
NMics = length(Vocs(1).Sound);
if NMics == 1; P.ComputeLocation = 0; end
% COMPUTE REFERENCE INDICES FOR SNR COMPUTATION
F = Vocs(1).F; RefInd = find(F>=20000 & F<=30000); 
% SELECT RANGE OF FREQUENCIES TO ANALYZE
FreqSelInd = find(F>=P.FreqRange(1) & F<=P.FreqRange(2));

%% FIND PHRASES
PhraseNumberByVoc = zeros(size(Vocs));
PhraseLengthByVoc = zeros(size(Vocs));
PhraseCounter = 1; PhraseLength = 0;
for iV=1:length(IVIs)
  if IVIs(iV) < P.WithinPhrasePause 
    PhraseLength = PhraseLength+1;
    PhraseNumberByVoc(iV) = PhraseCounter;
  else
    PhraseCounter = PhraseCounter + 1;
    PhraseNumberByVoc(iV) = PhraseCounter;
    PhraseLengthByVoc(iV-PhraseLength:iV-1) = PhraseLength;
    PhraseLength = 1;
  end
end

DeleteInd = [];

% Used for USM4 filtering before, also used in USVLoc2D Paper. Leads to better overall accuracy, but reduces localization at high mean Frequencies
%Filter = designfilt('bandpassiir','FilterOrder',20, 'HalfPowerFrequency1',30000,'HalfPowerFrequency2',100000, 'SampleRate',Vocs(1).SRSound);
%Filter = designfilt('highpassiir','FilterOrder',10, 'PassbandFrequency',35000,'SampleRate',Vocs(1).SRSound);
Filter = [];

%% PROCESS VOCALIZATIONS
for iV= P.Selection  % LOOP OVER VOCALIZATIONS
  if P.Verbose printupdate([num2str(iV),'/',num2str(length(P.Selection))],iV==1); end
   
  %% EXTRACT PROPERTIES FROM THE SPECTROGRAM
  SRSpec = Vocs(iV).SRSpec;
  
  % COMPUTE AVERAGE SPECTRAL ENERGY
  Vocs(iV).MeanEnergy = zeros(1,NMics); Vocs(iV).MaxEnergy = zeros(1,NMics);
  for iM =1:NMics
    Vocs(iV).MeanEnergy(iM) =  mean(mat2vec(Vocs(iV).Spec{iM}(FreqSelInd,:)));  
    Vocs(iV).MaxEnergy(iM) =  mean(max(Vocs(iV).Spec{iM}(FreqSelInd,:)));  
  end
  
  % SELECT SPECTROGRAM FOR FURTHER ANALYSIS WITH BEST SNR
  iLarger = NaN;
  switch NMics
    case 1;  iLarger = 1; cSpecOpt = Vocs(iV).Spec{iLarger};
    otherwise 
      % FIND MICROPHONE WITH LOUDEST MAXIMAL BIN
      [~,SortInd] = sort(Vocs(iV).MaxEnergy,'descend');
      iLarger = SortInd(1); cSpecOpt = Vocs(iV).Spec{iLarger(1)};
      if Vocs(iV).MaxEnergy(SortInd(2)) > 0.8*Vocs(iV).MaxEnergy(SortInd(1))
        iLarger(2) = SortInd(2); cSpecOpt = (cSpecOpt + Vocs(iV).Spec{iLarger(2)})/2;
      end
      
  end
  cSpecOpt(setdiff([1:end],FreqSelInd),:) = 0;
  %cSpecOpt(cSpecOpt<0.15 * max(cSpecOpt(:))) = 0;
  
  % FILTER SPECTRUM BEFORE FURTHER PROCESSING
  switch P.SpectrogramFilter
    case 'CurvatureFlow';
      cSpecOpt = curvatureFlowFilter(cSpecOpt,5,0.5,0.05,30);
      %cSpecOpt = MinMax_Noise_Removal(cSpecOpt,3,0.3,0.05,30); % original code, quite slow.
    case 'Wiener';
      cSpecOpt = wiener2(cSpecOpt,[5,5],0.5);
    case 'None'
    otherwise error('Method not implemented');
  end
  
  % COMPUTE TIME STEPS OF ACTUAL VOCALIZATION
  % SPECTRUM
  E = edge(cSpecOpt,'Prewitt',P.EdgeThreshold);
  E_Sum = sum(E);
  VocInds = find(E_Sum(P.BorderBins+1:end-P.BorderBins)>0);
  if length(VocInds)<3
    if P.ExcludeNonUSVs
      DeleteInd(end+1) = iV; continue;
    else
      VocInds = [1:3];
    end
  end
  VocRange = [VocInds(1) : VocInds(end)] + P.BorderBins;
    
  % REASSIGN TIME RANGE
  StartTime = VocRange(1)/SRSpec;
  StopTime = VocRange(end)/SRSpec;
  Vocs(iV).StartTime = StartTime;
  Vocs(iV).StopTime = StopTime;
  Vocs(iV).Duration = StopTime-StartTime;
  
  % SUBSELECT SPECTRUM OF VOCALIZATION
  cSpecOpt = cSpecOpt(:,VocRange);
  VocInds = VocInds - VocInds(1)+1;
  
  % COMPUTE THE SPECTRAL LINE (= frequency with maximal intensity in every bin)
  cSpecLine = NaN*zeros(1,length(VocRange));
  cSpecLineInd = NaN*zeros(1,length(VocRange));
  for iP = 1:length(VocInds) % STEP ALONG TIME
    iT = VocInds(iP);
    if any(cSpecOpt(:,iT))
      [~,cPos] = max(cSpecOpt(:,iT)); cSpecLine(iT) = F(cPos);
      cSpecLineInd(iT) = cPos;
    end
  end
  
  % CORRECT SPECTRAL LINE FOR NOISE
  NSteps = length(VocRange); Inds = [];
  for iP = [1,2,3:NSteps-2,NSteps,NSteps-1] % [floor(NSteps/2):-1:1,ceil(NSteps/2):NSteps]; % Work from inside to outside
    switch iP
      case {1,2};   %START & END
        cMean = cSpecLine([3]);
        if abs(1  - cSpecLine(iP)/cMean ) > 0.05
          cSpecLine(iP) = cMean;
        end
        continue;
      case {length(VocRange),length(VocRange)-1}
        cMean = cSpecLine([length(VocRange)-2]);
        if abs(1  - cSpecLine(iP)/cMean ) > 0.05
          cSpecLine(iP) = cMean;
        end        
        continue; 
      otherwise Inds = [-1,1];  % MIDDLE
    end
    % AVERAGE OVER NEIGHBORS IF DEVIATION TOO STRONG
    if length(cSpecLine) >= 3
      cMean = nanmean(cSpecLine(iP+Inds));
      IsolatedNaN = isnan(cSpecLine(iP)) & ~isnan(cSpecLine(iP-1)) & ~isnan(cSpecLine(iP+1)); % if isolated NaN between non-NaNs
      BigLocalChange = abs(1  - cSpecLine(iP)/cMean ) > 0.05; % if Change in one ms greater than 5%
      SimilarNeighbors= abs(1 - cSpecLine(iP+Inds(1))/cSpecLine(iP+Inds(2)))<0.1; % Less than 10% difference between neighboring bins
   
      if SimilarNeighbors & (IsolatedNaN | BigLocalChange)
        cSpecLine(iP) = cMean;
        [~,cSpecLineInd(iP)] = min(abs(F-cMean));
      end
    end
  end

  Vocs(iV).SpecLine = cSpecLine;
  
  %% COMPUTE BEST MICROPHONE BY TAKING THE SNR OF THE SPECTRAL LINE IN RELATION TO OTHER BINS
  VocInd = zeros(1,size(cSpecOpt,2));
  for iP=1:size(cSpecOpt,2)
    [~,VocInd(iP)] = min(abs(cSpecLine(iP)-F));
  end
  LinInd = sub2ind(size(cSpecOpt),VocInd,[1:size(cSpecOpt,2)]);
  Vocs(iV).SNRVoc = zeros(1,NMics);
  for iM = 1:NMics
    ccSpec = Vocs(iV).Spec{iM}(:,VocRange);
    cMean = nanmean(ccSpec(LinInd));
    if cMean > 0
      Vocs(iV).SNRVoc(iM) = cMean/nanmean(nanmean(ccSpec(RefInd,:)));
    else
      Vocs(iV).SNRVoc(iM) = 0;
    end
  end
  assert(sum(isnan(Vocs(iV).SNRVoc))==0,'SNR_Voc computation failed!');
  
  %% COMPUTE MANY BASIC PROPERTIES
  % RECOMPUTE SPECTRAL PURITY
  SD = NaN*zeros(size(cSpecOpt,2),1); 
  cSpecOptMax = max(cSpecOpt);
  for iT=1:size(cSpecOpt,2) % STEP THROUGH TIME
    cMed = median(cSpecOpt(:,iT));
    if cMed > 0; SD(iT) = cSpecOptMax(iT)/cMed; end
  end
  Vocs(iV).SpecPurity = nanmean(SD);
  % CENTER FREQUENCY
  Vocs(iV).FMean = nanmean(cSpecLine);
  % MAXIMAL FREQUENCY
  Vocs(iV).FMax = prctile(cSpecLine,90);
  % MINIMAL FREQUENCY
  Vocs(iV).FMin = prctile(cSpecLine,10);
  % FREQUENCY RANGE
  Vocs(iV).FRange = Vocs(iV).FMax-Vocs(iV).FMin;
  % STARTING FREQUENCY
  Vocs(iV).FStart = mean(cSpecLine(1:min(2,end)));
  % STOPPING FREQUENCY
  Vocs(iV).FStop = mean(cSpecLine(max(1,end-1):end));
  % DIRECTIONALITY
  Vocs(iV).Directionality = nanmean(diff(cSpecLine));
  % MARGINAL AMPLITUDE DISTRIBUTION OF VOCALIZATION in TIME
  Vocs(iV).Marginal.Time = NaN*zeros(size(cSpecLine));
  for iT = 1:length(cSpecLineInd)
    if ~isnan(cSpecLineInd(iT))
      Vocs(iV).Marginal.Time(iT) = cSpecOpt(cSpecLineInd(iT),iT);
    end
  end
  % SKEW AND KURTOSIS OF THE TIME MARGINAL
  cTime = [1:length(Vocs(iV).Marginal.Time)]/Vocs(iV).SRSpec;
  Vocs(iV).Marginal.TimeSkew = skewnesspdf(cTime,Vocs(iV).Marginal.Time);
  Vocs(iV).Marginal.TimeKurtosis = kurtosispdf(cTime,Vocs(iV).Marginal.Time);
  
  % MARGINAL AMPLITUDE DISTRIBUTION OF VOCALIZATION in FREQUENCY
  Vocs(iV).Marginal.Freq = zeros(size(F));
  for iT=1:length(cSpecLineInd)
    iF = cSpecLineInd(iT);
    if ~isnan(iF)
      Vocs(iV).Marginal.Freq(iF) = Vocs(iV).Marginal.Freq(iF) + Vocs(iV).Marginal.Time(iT);
    end
  end
  Vocs(iV).Marginal.Freq = Vocs(iV).Marginal.Freq/length(cSpecLine);
  % SKEW AND CURTOSIS OF THE FREQUENCY MARGINAL
  Vocs(iV).Marginal.FreqSkew = skewnesspdf(F,Vocs(iV).Marginal.Freq);
  Vocs(iV).Marginal.FreqKurtosis = kurtosispdf(F,Vocs(iV).Marginal.Freq);
  
  % COMPUTE THE WIENER ENTROPY
  % geometric mean, divided by arithemetic mean
  WienerEntropies = zeros(1,size(cSpecOpt,2));
  Tmp = cSpecOpt(2:end-1); Tmp(Tmp==0) = min(Tmp(Tmp>0))-eps;
  Vocs(iV).WienerEntropy = mean(geomean(Tmp)/mean(Tmp));
  
  % COMPUTE THE PITCH SALIENCE
  % https://essentia.upf.edu/reference/streaming_PitchSalience.html
  if P.PitchSalience
    cPS = zeros(NSteps,1); NFStepsMax = 50;
    for iT=1:NSteps
      cX = xcorr(D(:,iT),NFStepsMax,'unbiased');
      cX = cX(NFStepsMax+1:end);
      cPos = findLocalExtrema(cX,'min',1);
      cPS(iT) = max(cX(cPos:end))/cX(1);
    end
    Vocs(iV).PitchSalience = mean(cPS);
  end
  
  % ESTIMATE VIBRATO
  % Essentially PitchSalience of the SpectralLine
  NFStepsMax = 20;
  PSpecLine = fft(cSpecLine);
  cX = xcorr(abs(PSpecLine),NFStepsMax,'unbiased');
  cX = cX(NFStepsMax+1:end);
  cPos = findLocalExtrema(cX,'min',1);
  if isempty(cPos) cPos = length(cX); end
  Vocs(iV).Vibrato = nanmax(cX(cPos:end))/cX(1);
  
  % ESTIMATE AMOUNT OF LOW FREQUENCY NOISE
  Vocs(iV).Spectrum = zeros(size(Vocs(iV).Spec{1},1),1);
  for iM = 1:NMics
    Vocs(iV).Spectrum = Vocs(iV).Spectrum + max(Vocs(iV).Spec{iM},[],2);
  end
  
  % VARIANCE 
  Vocs(iV).Variance = var(cSpecLine);
  % VARIANCE OF SPECTRAL LINE
  Vocs(iV).SpecLineVar = var(diff(cSpecLine));
  % LOCAL VARIABILITY
  Vocs(iV).LocalChange = sum(abs(diff(cSpecLine)))/size(Vocs(iV).Spec,2);   
  
  % POSITION IN CURRENT PHRASE AND PHRASE NUMBER
  LocalPosition = find(IVIs(iV:-1:1)>0.2,1,'first');
  if isempty(LocalPosition) LocalPosition = 1; end
  Vocs(iV).LocalPosition = LocalPosition;
  Vocs(iV).PhraseNumber = PhraseNumberByVoc(iV);
  Vocs(iV).PhraseLength = PhraseLengthByVoc(iV);
 
  %% EXTRACT PROPERTIES FROM THE SOUND PRESSURE DATA
  SRSound = Vocs(iV).SRSound;
  SoundSteps = length(Vocs(iV).Sound{1});
  iSoundStart = round(StartTime*SRSound)+1;
  iSoundStop = min(round(StopTime*SRSound), length(Vocs(iV).Sound{1}));
  
  % COMPUTE AVERAGE LEVEL BASED ON SOUND
  for iM=1:NMics
    Vocs(iV).Level(iM) = std(Vocs(iV).Sound{iM}(iSoundStart:iSoundStop));
    Vocs(iV).BaselineStd(iM) = std(Vocs(iV).Sound{iM}([1:iSoundStart-1,iSoundStop:end]));
    Vocs(iV).SNR(iM) = Vocs(iV).Level(iM)/Vocs(iV).BaselineStd(iM);
  end
  
  Vocs(iV).BaselineStd(iM) = std(Vocs(iV).Sound{iM}(1:iSoundStart));
  Vocs(iV).SNR(iM) = Vocs(iV).Level(iM)/Vocs(iV).BaselineStd(iM); 
  
  %% ESTIMATE SPATIAL LOCATION OF THE VOCALIZATION FROM MULTIPLE MICROPHONES
  if P.ComputeLocation
    StartTime = Vocs(iV).StartTime;
    StopTime = Vocs(iV).StopTime;
    SR = Vocs(iV).SRSound;
    LocationBinning = 'MultiWindow';
    VSound = 343; % [m/s]
    CorrTime = P.CorrRange/VSound;
    TimeSteps = []; Range = [];
    switch LocationBinning
      case 'SingleWindow'
        TimeSteps = Vocs(iV).Duration/2;
        Range = max([Vocs(iV).Duration/2,0.03]);
      case 'MultiWindow'
        TimeSteps =  [0:0.003:Vocs(iV).Duration];
        Range = 0.03; % s, before and after to analyze
    end
    % Modify to only start at the beginning of the vocalization and 
    LocTmp = zeros(length(TimeSteps),2);
    CertTmp = LocTmp; S2NTmp = LocTmp(:,1);
    for iPos=1:length(TimeSteps)
      RelTime  = TimeSteps(iPos);
      StartStep = round(SR * (StartTime + RelTime - Range));
      StopStep = round(SR * (StartTime + RelTime + Range));
      NSteps    = length(Vocs(iV).Sound{iM});
      Ind           = [max([1,StartStep]) : min([NSteps,StopStep]) ];
      LocSound = cell(NMics,1);
      for iM = 1:NMics
        LocSound{iM} = Vocs(iV).Sound{iM}(Ind);
      end
      RLoc = ...
        VocLocalizer('Sounds',LocSound,'SR',Vocs(iV).SRSound,'Filter',Filter,...
        'CorrMethod',{'EWGCC'},'CenterShift',P.CenterShift,...
        'MicrophonePositions',P.MicrophonePositions,'SourceHeight',P.SourceHeight,...
        'EstimationMethod',P.EstimationMethod,'FIG',0,'CorrTime',[]);
      LocTmp(iPos,1:2) = RLoc.StimPosMean;
      CertTmp(iPos,1:2) = RLoc.StimPosSD;
      S2NTmp(iPos,1) = nanmean(RLoc.SignalToNoiseAll(:));
    end
    % FOR DIMENSIONS > 1, THE CERTAINTY CAN BE USED
    CertTmp = max(CertTmp,[ ],2);
    cInd = find(CertTmp<0.01); % Only keep certainties below a centimeter
    if ~isempty(cInd)
      Vocs(iV).Location = median(LocTmp(cInd,:),1);
      Vocs(iV).LocationCertainty = median(CertTmp(cInd),1);
    else 
      Vocs(iV).Location = median(LocTmp,1);
      Vocs(iV).LocationCertainty = median(CertTmp,1);
    end
    % FOR ALL DIMENSIONS, THE SIGNAL TO NOISE OF THE XC CAN BE USED
    Vocs(iV).SignalToNoise = median(S2NTmp);
    
    if sum(imag(Vocs(iV).Location)) keyboard; end
  end
  
  %% PLOT PROPERTIES FOR CHECKING
  if mod(iV,P.CheckProperties)==0; LF_showVoc(Vocs(iV),NMics,VocRange,cSpecOpt); end
  
end
fprintf('\n')

Vocs = Vocs(setdiff(1:length(Vocs),DeleteInd));

%% SELECT BASED ON CERTAIN PROPERTIES
IndAll = logical(ones(size(Vocs)));
for iM = 1:length(P.Selector)
  Property = P.Selector{iM}{1};
  Condition = P.Selector{iM}{2};
  Values = [Vocs.(Property)];
  cInd = eval(['Values ',Condition]);
  IndAll = logical(IndAll.*cInd);
end
Vocs = Vocs(IndAll);


%% SHOW VOCALIZATIONS
function LF_showVoc(Voc,NMics,VocRange,SpecOpt)

% PREPARE FIGURE
figure(1001); clf;
set(gcf,'PaperSize',[8.5,10],'PaperPosition',[0,0,8.5,10]);
AxisLabelOpt = {'FontSize',8};
[DC,AH] = axesDivide(ones(1,NMics+1),[2,1],[0.16,0.12,0.8,0.8],'c');
for iA=1:numel(AH) set(AH(iA),'FontSize',7); end

% PREPARE PLOTTING
TimeSound = ([1:length(Voc.Sound{1})]/Voc.SRSound - Voc.StartTime)*1000;
TimeSpec = ([1:size(Voc.Spec{1},2)]/Voc.SRSpec - Voc.StartTime)*1000;
CMAX = max(mat2vec(cell2mat(Voc.Spec)));
for iM=1:NMics+1
  % VOLTAGE
  if iM<=NMics
    axes(AH(2,iM)); hold on;
    PreSteps = [1:round(Voc.StartTime * Voc.SRSound)];
    SD = std(Voc.Sound{iM}(PreSteps));
    plot(TimeSound,Voc.Sound{iM}/SD,'k');
    xlim([TimeSound([1,end])]);
    xlabel('Time [ms]',AxisLabelOpt{:});
    ylim([-10,10]);
    if iM==1; ylabel('Voltage [noise S.D.]',AxisLabelOpt{:}); end
  end
  
  % SPECTROGRAM
  axes(AH(1,iM));  hold on; box on;
  if iM<=NMics cSpec = Voc.Spec{iM}; else cSpec = SpecOpt; end
  imagesc(TimeSpec,Voc.F/1000,cSpec);
  set(gca,'YDir','normal');
  xlim([TimeSpec([1,end])]+[-5,5]);
  ylim([0,125]);
  Range = TimeSpec([1,end]);
  plot(Range,repmat(Voc.FMin/1000,[2,1]),'b');
  plot(Range,repmat(Voc.FMean/1000,[2,1]),'Color',[0,0.5,1],'LineWidth',2);
  plot(Range,repmat(Voc.FMax/1000,[2,1]),'b');
  plot(0,Voc.FStart/1000,'.','Color',[0,0,1],'MarkerSize',14);
  plot(Voc.Duration*1000,Voc.FStop/1000,'.','Color',[0.2,0.2,1],'MarkerSize',14);
  plot(TimeSpec(VocRange),Voc.SpecLine/1000,'r','LineWidth',2);
  if iM == 1; ylabel('Freq. [kHz]',AxisLabelOpt{:});  end
  if iM <= NMics
    title(['Microphone ',num2str(iM),' (SNRVoc = ',num2str(Voc.SNRVoc(iM),2),')'],'FontSize',8,'FontWeight','bold');
  end
  colormap(HF_colormap({[1,1,1],[0,0,0]},[0,1]));
  if iM ==1 text(0.02,0.04,[Voc.Animal,'r',num2str(Voc.Recording),'t',num2str(Voc.Start,4)],'FontSize',5,'Units','n'); end
  caxis([0,CMAX/1.5]);
end
pause