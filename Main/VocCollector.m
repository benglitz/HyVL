function Vocs = VocCollector(varargin)
% 
% Example Usage:
%  Vocs = VocCollector('Animals',{'mouse3'},'Recording',80);
% 
%  Vocs = VocAnalyzer(Vocs);
  
P = parsePairs(varargin);
checkField(P,'DataSource','Controller'); % CAN BE EITHER Controller, or WAV, or Direct
checkField(P,'Data',''); % If DataSource = Direct, then the Data should be passed here.
% CONTROLLER PARAMETERS
checkField(P,'Paradigm','');
checkField(P,'Animals',{});
checkField(P,'Recording',[]);
% WAV PARAMETERS
checkField(P,'Filename',[]);
if strcmp(P.DataSource,'WAV') & isempty(P.Filename);
  error('Filename has to be provided if DataSource is WAV'); end;
% EXTRACTION PARAMETERS
checkField(P,'VocJitter',0.002);
checkField(P,'FWindow',0.001);
checkField(P,'FRange',[5000,125000]);
checkField(P,'Highpass',0);
checkField(P,'PreTime',0);
checkField(P,'PostTime',0);
checkField(P,'Channels',[1:4]);
checkField(P,'Reload',1);
checkField(P,'SRAI',[]);
checkField(P,'Classifications',[]);
checkField(P,'Plot',0);
checkField(P,'TVocs',[]); % Externally determined times of vocalizations
checkField(P);

if ~iscell(P.Animals) P.Animals = {P.Animals}; end
if ~isempty(P.Data) P.DataSource = 'Direct'; end

switch P.DataSource
  case 'Controller';
    % SELECT RECORDINGS FROM DATA BASE
    R = [];
    %C_checkDatabase('Force',1);
    % GET ANIMAL IDs
    RA = mysql('SELECT * FROM subjects');
    for i=1:length(RA) IDsByAnimal.(RA(i).name) =RA(i).id; end
    if isempty(P.Animals) P.Animals = setdiff({RA.name},'test'); end
    
    % GET RECORDING INFO
    for iA=1:length(P.Animals)
      MYSQL = ['SELECT * FROM recordings WHERE '...
        'animal=''',P.Animals{iA},''' AND bad=0  '];
      if ~isempty(P.Recording)
        MYSQL = [MYSQL,' AND recording=',num2str(P.Recording)];
      end
      if ~isempty(P.Paradigm)
        MYSQL = [MYSQL,' AND paradigm=''',P.Paradigm,''''];
      end
      R = [R;mysql(MYSQL)];
    end
    
    fprintf(['\n= = = = Found [ ',num2str(length(R)),' ] Recordings = = = =\n\n']);
    % MAKE RECORDINGS LOCAL
    for iR=1:length(R)
      Path = C_makeLocal('Animal',R(iR).animal,'Recording',R(iR).recording,...
        'Modules',{'AnalogIn','NIDAQ'});
    end
    
  case 'Local'
    R = struct('animal',P.Animals,'recording',P.Recording);
    
  case 'WAV' % DEFINE RECORDINGS
    R.Filename = P.Filename; % Just one filename possible
    
  case 'Direct' % 
    R.Filename = 1;
    
end
    
k=0; Vocs = [];
for iR=1:length(R)
  
  % LOAD DATA VOCALIZATIONS
  global CurrentData
  if isempty(CurrentData) || P.Reload || length(R)>1
    switch P.DataSource
      case {'Controller','Local'}; 
        fprintf(['\n= = = = = =  Recording : [ ',R(iR).animal,' R',num2str(R(iR).recording),' ] = = = = = = = = = = = =\n']);
        CurrentData = C_loadRecording('Animal',R(iR).animal,'Recording',R(iR).recording,...
          'Modules',{'AnalogIn','NIDAQ'});
        SRAI = CurrentData.General.Parameters.Setup.Audio.SRAI;

      case 'WAV';
        fprintf(['\n= = = = = =  Recording : [ ',escapeMasker(P.Filename),' ] = = = = = = = = = = = =\n']);
        [Data,SRAI] = audioread(P.Filename);
        CurrentData.AnalogIn.Data(1).Data.Analog = Data;
        CurrentData.AnalogIn.Data(1).Data.Time = 0;
        P.Channels = intersect(P.Channels,[1:size(Data,2)]);
        
      case 'Direct'
        CurrentData.AnalogIn.Data(1).Data.Analog = P.Data;
        CurrentData.AnalogIn.Data(1).Data.Time = 0;
        P.Channels = [1:size(P.Data,2)];
        SRAI = P.SRAI;
        
      otherwise error('Data Source not implemented.');
        
    end
    
    Reloaded = 1;
  else Reloaded = 0;
  end
  
  CurrentData.AnalogIn.Data(1).Data.Analog = CurrentData.AnalogIn.Data(1).Data.Analog(:,P.Channels);
  NChannels = size(CurrentData.AnalogIn.Data(1).Data.Analog,2);

  % DOWN SAMPLE FOR TOO HIGH SAMPLE RATES
  if SRAI == 500000 && Reloaded
    CurrentData.AnalogIn.Data(1).Data.Analog =  CurrentData.AnalogIn.Data(1).Data.Analog(1:2:end,:);
    CurrentData.AnalogIn.Data(1).Data.Time =  CurrentData.AnalogIn.Data(1).Data.Time(1:2:end);
    SRAI = 250000;
  end
  
  % HIGHPASS FILTER TO DEEMPHASIZE ENVIRONMENTAL NOISE
  if P.Highpass
    [b,a] = butter(4,[P.Highpass]/(SRAI/2),'high');
    warning('filtfilt not tested yet. Preferable since this allows zero-phase filtering');
    CurrentData.AnalogIn.Data(1).Data.Analog =  filtfilt(b,a,CurrentData.AnalogIn.Data(1).Data.Analog);
  end 
  
  % EXTRACT VOCALIZATIONS
  RecordingOnset = CurrentData.AnalogIn.Data(1).Data.Time(1);
  FWindowSteps = P.FWindow*SRAI;
  FWindowSteps = 500;
  if NChannels > 0 % DATA NOT MISSING
    for iT=1:length(CurrentData.AnalogIn.Data) % LOOP OVER TRIALS
      fprintf(['    Trial ',num2str(iT),'\n']);
      
      for iM=1:NChannels % LOOP OVER MICROPHONES
        fprintf(['    > Computing Spectrogram on Channel ',num2str(iM),'\n']);
        % COMPUTE SPECTROGRAM
        [S{iM},F,T,Thresh(iM),SRSpec] = HF_specgram(CurrentData.AnalogIn.Data(iT).Data.Analog(:,iM),...
          FWindowSteps,SRAI,P.FRange,FWindowSteps/2,1,1);
        
        % EXTRACT VOCALIZATIONS
        if isempty(P.TVocs);
          TVocs{iM} = LF_findVocs(S{iM},F,T,'Threshold',Thresh(iM),'Channel',iM,'Plot',P.Plot);
        else % Use externally provided ones if available
          TVocs{iM} = P.TVocs;
        end
      end
      
      if ~isempty(cell2mat(TVocs))
        if NChannels == 1; TVocsAll = TVocs{1}; InterVocTimes = [NaN,diff(TVocsAll(1,:))];
        else  [TVocsAll,InterVocTimes] = LF_fuseVocTimes(TVocs,P.VocJitter,SRSpec);
        end
        
        [cVocsSound,cVocsSpec,TVocsAll] = LF_getVocs(iT,S,TVocsAll,SRAI,SRSpec,P.PreTime,P.PostTime);
        
        if ~isempty(TVocsAll)
          for iV=1:size(TVocsAll,2) % LOOP OVER FOUND VOCALIZATIONS
            k=k+1;
            for iM=1:NChannels % LOOP OVER MICROPHONES
              % COLLECT IN AN INTERMEDIATE FORMAT
              Vocs(k).Sound{iM} = cVocsSound{iM}{iV};
              Vocs(k).Spec{iM} = full(abs(cVocsSpec{iM}{iV}));
              Vocs(k).SpecPurity{iM} = LF_SpecPurity(Vocs(k).Spec{iM});
              switch P.DataSource
                case {'Controller','Local'};
                  Vocs(k).Recording = R(iR).recording;
                  Vocs(k).Animal = R(iR).animal;
                  Vocs(k).AnimalNumber = str2num(R(iR).animal(6:end));
                  switch P.Paradigm
                    case 'Interaction';
                      if isfield( CurrentData.General.Paradigm.Parameters,'InteractionPartner')
                        Vocs(k).AnimalPartner = CurrentData.General.Paradigm.Parameters.InteractionPartner;
                        Vocs(k).AnimalPartnerNumber = str2num(Vocs(k).AnimalPartner);
                      end
                  end
                case 'WAV';
                  Vocs(k).Filename = R(iR).Filename;
                  Vocs(k).Animal = 'Unknown'; Vocs(k).Recording = 0;
              end         
              Vocs(k).Trial = iT;
              Vocs(k).Number = iV;
              % START AND STOP OF CUTOUT DATA
              Vocs(k).StartWin = TVocsAll(1,iV) + RecordingOnset;
              Vocs(k).StopWin = TVocsAll(2,iV) + RecordingOnset;
              Vocs(k).DurationWin = Vocs(k).StopWin - Vocs(k).StartWin;
              % START AND STOP OF ACTUAL VOCALIZATION
              Vocs(k).Start = TVocsAll(1,iV)+P.PreTime + RecordingOnset;
              Vocs(k).Stop = TVocsAll(2,iV)- P.PostTime + RecordingOnset;
              Vocs(k).Duration = Vocs(k).Stop - Vocs(k).Start;
              
              Vocs(k).PreTime = P.PreTime;
              Vocs(k).PostTime = P.PostTime;
              Vocs(k).Interval = InterVocTimes(iV);
              Vocs(k).SRSound = SRAI;
              Vocs(k).SRSpec = SRSpec;
              Vocs(k).dF = F(2)-F(1);
              Vocs(k).F = F;
              for iC=1:length(P.Classifications)/2
                Vocs(k).(P.Classifications{(i-1)*2+1}) = P.Classifications{i*2};
              end
              Vocs(k).Time = Vocs(k).StartWin + [1:size(Vocs(k).Spec{iM},2)]/Vocs(k).SRSpec;
            end
          end
        end
      end
    end
  end
  if ~isempty(Vocs)
    fprintf(['Minimal Duration: ',num2str(min([Vocs.DurationWin])),'\n']);
    fprintf(['Average Duration: ',num2str(mean([Vocs.DurationWin])),'\n']);
    fprintf(['Maximal Duration: ',num2str(max([Vocs.DurationWin])),'\n']);
  end
end
clear global CurrentData
  
% FIND THE VOCALIZATIONS
function TVocs = LF_findVocs(S,F,T,varargin)
    
    P = parsePairs(varargin);
    checkField(P,'Channel');
    checkField(P,'DurationMin',0.007); % Vocalizations have to be longer than 7ms
    checkField(P,'DurationMax',0.500); % Vocalizations have to be longer than 
    checkField(P,'FreqMeanThresh',35000);
    checkField(P,'FreqLowThresh',1.5); % Check that the vocalization does not have a lot of energy in very low frequencies
    checkField(P,'MaxEnergyThresh',1); % This Threshold is not in relative units!
    checkField(P,'SpecPurityThresh',0);
    checkField(P,'SpecContThresh',20);
    checkField(P,'SpecDiscThresh',0); % 0.7?
    checkField(P,'MergeClose',0.015);
    checkField(P,'FilterDuration',0.01);
    checkField(P,'Plot',0);
    
    fprintf(['\tAnalyzing Channel ',num2str(P.Channel),'  : ']);
   
    dT = diff(T(1:2)); dF = diff(F(1:2));
    FilterSteps = round(P.FilterDuration/dT);
    [NFreq,NTime] = size(S);
    iBad = zeros(0,NTime); Criteria = {};
    
    % CONVERT SPECGRAM TO SPECTRAL POWER  
    S = sparse(medfilt2(full(S),[1,5]));
    [i,j,s] = find(S);   SPow = S.^2;
    totPower = sum(SPow);   
    iNonZero = find(totPower);
    
    % POWER THRESHOLD 
    PowerBins = [0:0.01:100];
    H = histc(full(SPow(:)),PowerBins);
    PowerThreshold = PowerBins(find(cumsum(H)/sum(H)>0.998,1,'first'));
    PowerThreshold = max([PowerThreshold,0.1]);
    
    % SELECT BY MAXIMAL POWER THRESHOLD
    if P.MaxEnergyThresh
      maxEnergy = max(S);
      iBad(end+1,:) = maxEnergy < P.MaxEnergyThresh;
      Criteria{end+1} = 'MaxEnergy';
      Values.(Criteria{end}) = maxEnergy;
    end
      
    % SELECT BY SPECTRAL CONTINUITY
    % Spectral Continuity : lower values are more continuous
    if ~isempty(P.SpecContThresh)
      SpecCont = C_SpectralContinuity(SPow,PowerThreshold);
      iBad(end+1,:) = SpecCont > P.SpecContThresh;
      Criteria{end+1} = 'SpecCont';
      Values.(Criteria{end}) = SpecCont;
    end
    
    % SELECT BY SPECTRAL DISCONTINUITY
    % Spectral Discontinuity
    if P.SpecDiscThresh
      SD = specdiscont(SPow);
      SD = medfilt1(SD,FilterSteps);
      iBad(end+1,:) = SD > P.SpecDiscThresh;
      Criteria{end+1} = 'SpecDisc';
      Values.(Criteria{end}) = SD;
    end

    % SELECT BY LOW FREQUENCY CONTENT
    if P.FreqLowThresh
      cIndLow = find(F<25000);
      cIndHigh = find(F>25000);
      FreqLow = mean(S(cIndLow,:))./mean(S(cIndHigh,:));
      iBad(end+1,:) = FreqLow > P.FreqLowThresh;
      Criteria{end+1} = 'FreqLow';
      Values.(Criteria{end}) = FreqLow;
    end
    
    % SELECT BY AVERAGE FREQUENCY 
    % mostly excludes loud low frequency noise
    if P.FreqMeanThresh;
      spfreq = sparse(i,j,i,size(S,1),size(S,2));
      MeanFreq = sum(spfreq.*SPow)*dF;
      MeanFreq(iNonZero) = MeanFreq(iNonZero)./totPower(iNonZero);
      if P.FilterDuration
        MeanFreq = medfilt1(full(MeanFreq),FilterSteps);
      end
      iBad(end+1,:) = MeanFreq < P.FreqMeanThresh;
      Criteria{end+1} = 'FreqMean';
      Values.(Criteria{end}) = MeanFreq;
    end
    
    % SELECT BY SPECTRAL PURITY
    % spectral purity evaluates how localized the sound is in a given frequency bin
    % Currently not used, since this is problematic for vocalization with wide, complex spectrogram
    if P.SpecPurityThresh  
      SpecPurity = zeros(size(T));
      SpecPurity(iNonZero) = maxPower(iNonZero)./totPower(iNonZero);
      
      if isequal(P.PurityThresh,'auto')
        [H,B] = hist(SpecPurity,[100]);
        [~,Ind] = max(H); BasePurity = B(Ind);
        Ind = find(SpecPurity<BasePurity);
        SD = sqrt(sum(SpecPurity(Ind) - BasePurity).^2)/length(Ind);
        P.PurityThresh = BasePurity + 5*SD;
      end
      
      if P.FilterDuration
        SpecPurity = medfilt1(SpecPurity,FilterSteps);
      end
      iBad(end+1,:) = SpecPurity <= P.SpecPurityThresh;
      Criteria{end+1} = 'SpecPurity';
      Values.(Criteria{end}) = SpecPurity;
    end
    
    % FIND INDICES THAT FULLFILL ALL CRITERIA
    iBadInd = find(max(iBad,[],1)); % transform to indices
    DurationMinSteps = P.DurationMin/dT;      iMin = find(diff(iBadInd) > DurationMinSteps);
    DurationMaxSteps = P.DurationMax/dT;     iMax = find(diff(iBadInd) < DurationMaxSteps);   
    iLength = intersect(iMin,iMax);
    TVocs = [T(iBadInd(iLength)); T(iBadInd(iLength+1))];
    fprintf([' ',num2str(size(TVocs,2)),' ==(merge close)==>']);
    
    % CHECK RESULTS
    if P.Plot
      figure(1); clf;
      set(gcf,'ButtonDownFcn','set(gca,''XLim'',5+get(gca,''XLim''))');
      TInd = 1:length(T);
      [~,AH] = axesDivide(1,[3,ones(1,length(Criteria))],'c');
      linkaxes(AH,'x');
      cT = T(TInd);
      axes(AH(1)); hold on
      imagesc(cT,F,S(:,TInd)); set(gca,'YDir','normal'); colorbar('Location','East');
      NVocs = size(TVocs,2);
      plot3(TVocs(1,:),repmat(80000,1,NVocs),repmat(1,1,NVocs),'.g'); 
      plot3(TVocs(2,:),repmat(80000,1,NVocs),repmat(1,1,NVocs),'.r'); 
      caxis([0,2]);
      for iC=1:length(Criteria)
        axes(AH(iC+1)); hold on;
        plot(cT,Values.(Criteria{iC})(TInd));
        cThresh = P.([Criteria{iC},'Thresh']);
        plot(cT([1,end]),[cThresh,cThresh],'r');
        ylabel(Criteria{iC})
      end
      set(AH,'XLim',[20,25]);
    end
    
    % MERGE CLOSE PAIRS
    if P.MergeClose
      DeltaT = TVocs(1,2:end)-TVocs(2,1:end-1);
      iClose = find(DeltaT <  P.MergeClose);
      for i = length(iClose):-1:1
        TVocs(2,iClose(i)) = TVocs(2,iClose(i)+1);
        TVocs(:,iClose(i)+1) = [];
      end
    end
    
    fprintf([' ',num2str(size(TVocs,2)),'\n']);
   
    
% FUSE THE VOCALIZATION TIMES ACROSS THE CHANNELS
function [TVocsAll,IVI] = LF_fuseVocTimes(TVocs,VocJitter,SRSpec)
  % FUSE THE WHISTLES ACROSS THE TWO SIDES
  fprintf('\tFusing Vocalizations  :  ');
  
  % CONSTRUCT COMMON VECTOR BETWEEN
  NChannels = length(TVocs);
  MaxDur = 0; 
  for iC=1:NChannels 
    if ~isempty(TVocs{iC})
      MaxDur = max(MaxDur,TVocs{iC}(2,end)); 
    end
  end
  MaxSteps = ceil(MaxDur*SRSpec);
  iGood = zeros(NChannels,MaxSteps);
  for iM=1:NChannels
    for iV = 1:size(TVocs{iM},2)
      cInd = [round(TVocs{iM}(1,iV)*SRSpec) : round(TVocs{iM}(2,iV)*SRSpec)];
      cInd = cInd(cInd>0);
      iGood(iM,cInd) = 1;
    end
  end
  
  % EXTRACT THE OVERLAPPING VOCALIZATIONS
  iFuse = [-1,find(sum(iGood,1)),MaxSteps+2];
  StopInd = find(diff(iFuse)>1);
  StopTimes = iFuse(StopInd(2:end))/SRSpec;
  StartTimes = iFuse(StopInd(1:end-1)+1)/SRSpec;
  TVocsAll = [StartTimes;StopTimes];
  IVI = [NaN,diff(TVocsAll(1,:))];
  NVocsAll = size(TVocsAll,2);
  fprintf([num2str(NVocsAll),'\n']);
  
% EXTRACT THE VOCALIZATIONS BASED ON THE TIMES
function [VSound,VSpec,TVocsAll] = LF_getVocs(iT,S,TVocsAll,SRSound,SRSpec,PreTime,PostTime)

  global CurrentData;
  
  NChannels = length(S);
  
  for iV = 1:size(TVocsAll,2)
    cTimes = TVocsAll(:,iV);
    cInd = round(cTimes*SRSpec);
    Started = 0;
    while ~Started
      for iC = 1:NChannels;   Started = abs(Started + S{iC}(:,cInd(1)));    end
      if ~Started cInd(1) = cInd(1)+1; end
    end
    Stopped = 0;
    while ~Stopped
      for iC = 1:NChannels;   Stopped = abs(Stopped + S{iC}(:,cInd(1)));    end
      if ~Stopped cInd(2) = cInd(2)-1; end
    end
    cTimes = cInd/SRSpec + [-PreTime;PostTime];
    TVocsAll(:,iV) = cTimes;
    cIndSpec = max(1,round(cTimes*SRSpec));
    cIndSound = max(1,round(cTimes*SRSound));
    for iM=1:NChannels
      % COLLECT SPECTROGRAM
      VSpec{iM}{iV} = S{iM}(:,cIndSpec(1):min(end,cIndSpec(2)));
      
      % COLLECT SOUND PRESSURE
      VSound{iM}{iV} = double(CurrentData.AnalogIn.Data(iT).Data.Analog(cIndSound(1):min(end,cIndSound(2)),iM));      
    end
  end
  
  function SpecPurity = LF_SpecPurity(SPow)
      totPower = sum(SPow);  
   [maxPower,maxIndx] = max(SPow);
    iNonZero = find(totPower);
    SpecPurity = zeros([size(SPow,2),1]);
    SpecPurity(iNonZero) = maxPower(iNonZero)./totPower(iNonZero);
    
 
