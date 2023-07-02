function MultiViewer(varargin)
% Data Viewer for Multimodal Data Collected with Controller
% Note: At this point MultiViewer only works with consecutive trials!
% 
% Usage: 
%  MultiViever('Argument',Value,...);
% For a list of arguments, see the code (edit Multiviewer) and see the list
% of checkField options at the beginning.
% 
% Keyboard Shortcuts :
% 
% Playing: 
% - Arrow left :  one frame back
% - Arrow right :  one frame forward
%
% Visual: 
% - g : show grid
% - +/- : change color range on the spectrogram 
% - S : Switch microphone
%
% TODO:
%
% DLC Tracking:
%   MultiViewer(39,43,'PhysicalSetup','Platform2D','SelTypes',{'Snout','LeftEar','RightEar','HeadCenter','TailStart'},'RandomFrames',100)
%
% Author : benglitz@gmail.com

%% PROCESS ARGUMENTS
if ~isempty(varargin) varargin = C_parseID(varargin); end
P =  parsePairs(varargin);
checkField(P,'Animal',[]);
checkField(P,'Recording',[]);
checkField(P,'Data',[]); % Multimodal data in Controller format (obtained from C_loadRecording)
% GENERAL OPTIONS
checkField(P,'PhysicalSetup',[]); % Choose the parameters of the setup the data was recorded with
checkField(P,'Lens','35mm');
checkField(P,'Time',[]); % Set the initially shown time
% VIDEO OPTIONS
checkField(P,'Transpose',0); % Transpose video (Rotate by 90 degree)
checkField(P,'Mirror1',0); % Mirror video horizontally
checkField(P,'Mirror2',0); % Mirror video vertically
checkField(P,'StepSize',0.1); % Size of time step for GUI arrows
checkField(P,'Window',2); % Time around the current time, displayed for Sensors, Ephys and Audio
checkField(P,'CMin',0); % Minimum of color range
checkField(P,'CMax',255); % Maximum of color range
checkField(P,'TrueLimits',0); % Compute true color range limits (slower startup)
checkField(P,'Annotations',[]);  % Annotations to display for each frame (for adding annotations)
checkField(P,'Lens','35mm'); % Lens used for the video recording
checkField(P,'Markers',{'o','^','v','d','x'}); % Marker shown for the different Tracking Points
checkField(P,'ShowGrid',1); %  Choose whether to show the grid in the video frame
checkField(P,'VideoShift',[0,0]);
% AUDIO OPTIONS
checkField(P,'SDAudio',NaN); % Set the standard deviation of the audio recording instead of estimating it (for speed) 
checkField(P,'ComputeLocation',0);
checkField(P,'Vocs',[]); % Use Vocalizations collected externally
checkField(P,'SelTypes',{'Snout','HeadCenter'});
checkField(P,'ShowZeroLine',1);
%checkField(P,'SelTypes',{'Snout','LeftEar','RightEar','HeadCenter','TailStart','Snout','LeftEar','RightEar','HeadCenter','TailStart'});
% Vocalizations :
checkField(P,'LoadAnnotations',1) % 
checkField(P,'SelTypesShift',{'Ascending','Descending','Ushaped','InvUshaped','Multipeaked','Break','Complex','Broadband','Other','Flat'});
checkField(P,'DetectVocs',0); % Automatically Detect Vocalizations
checkField(P,'StepTo','Annotation') % Options : Annotation / Vocalization
checkField(P,'NAudioChannels',[]);
checkField(P,'RandomFrames',[]);

% EPHYS OPTIONS
checkField(P,'Electrodes',[]); % Electrodes to show
checkField(P,'PerProngRef',1); % PerProngRefecence
checkField(P,'Reference','all');

% DISPLAY OPTIONS
checkField(P,'ShowSensors',1);   % Toggle Display of Sensors
checkField(P,'ShowAudio',1);     % Toggle Display of Audio
checkField(P,'ShowAnalogOut',1); % Toggle Display of AudioOut (not saved)
checkField(P,'ShowStimulus',1);  % Toggle Display of Audio
checkField(P,'ShowSpec',1);      % Toggle Display of Spectrogram
checkField(P,'ShowVideo',1);     % Toggle Display of Video
checkField(P,'ShowEPhys',0);     % Toggle Display of EPhys
checkField(P,'ShowSpikes',0);    % Toggle Display of Spikes
checkField(P,'CorrectDistortion',0); % Correct Video Distortion due to lens 
checkField(P,'FIG',1);           % Number of figure to plot into. 
checkField(P);

if ~isempty(P.Animal) && ~isempty(P.Recording)
  Modules = {'NIDAQ'};
  if P.ShowVideo; Modules{end+1} = 'Video'; end % Generalize over different Hardware: PointGrey or IDS
  if P.ShowAudio | P.ShowSpec; Modules{end+1} = 'AnalogIn'; end
  if P.ShowAnalogOut; Modules{end+1} = 'AnalogOut'; end
  if P.ShowEPhys; Modules{end+1} = 'EPhys'; end
  clear global R; R = C_loadRecording('Animal',P.Animal,'Recording',P.Recording, ...
      'VideoLoadMode','fullbuiltin','Modules',Modules,'Electrodes',P.Electrodes, ...
      'Reference', P.Reference, 'ShowStimulus', P.ShowStimulus); global R;
else
  warning off 
  if ~isempty(P.Data)
    global R; R = P.Data; rmfield(P,'Data');
  else
    if evalin('base','exist(''R'',''var'')')
      evalin('base','global R');
      global R;
    else
      fprintf('If Data is not passed via field ''Data'', a variable R needs to exist in the base workspace.\n'); return;
    end
  end
  warning on
  P.Animal = R.General.Parameters.General.Animal;
  P.Recording = R.General.Parameters.General.Recording;
end
checkField(P);

try delete(P.FIG); end
warning off; global MV R; MV = [];  warning on;

% CHOOSE VIDEO SOURCE
if isfield(R,'Thorcam') MV.VideoName = 'Thorcam';
elseif isfield(R,'PointGrey')   MV.VideoName = 'PointGrey'; 
elseif isfield(R,'IDSuEye') MV.VideoName = 'IDSuEye'; 
elseif isfield(R,'VideoCalcium') MV.VideoName = 'VideoCalcium'; 
elseif isfield(R,'Scanimage') MV.VideoName = 'Scanimage'; P.ShowGrid = 0;
else MV.VideoName = 'Video'; end

MV.P = P;
MV.Modules = setdiff( fieldnames(R),'General');

% GUESS CORRECT SETUP
if isempty(P.PhysicalSetup)
  cAnimal = R.General.Parameters.General.Animal;
  if strcmp(MV.VideoName,'Thorcam')
    P.PhysicalSetup = 'Widefield';
  else
    switch lower(cAnimal)
    case 'test'; fprintf('For animal ''test'' setup cannot be guessed.'); return;
    otherwise
      cNum = str2num(cAnimal(6:end));
      if any(cNum == [49:94]) ; P.PhysicalSetup = 'Platform2DMic4Precise'; end
      if any(cNum >= 95) ; P.PhysicalSetup = 'AudioBehavior'; end
      if isempty(P.PhysicalSetup) error('Could not guess PhysicalSetup. Provide during call of MultiViewer.'); end
    end
  end
end

if isstr(P.PhysicalSetup)
  MV.SetupInfo = C_getSetupInfo('PhysicalSetup',P.PhysicalSetup,'Lens',P.Lens,'VideoShift',P.VideoShift);
else 
  MV.SetupInfo = P.PhysicalSetup;
end
MV.FIG=P.FIG;
MV.Colors  = struct('NIDAQ',[0,0,0],'AnalogIn',[0,0,1],MV.VideoName,[1,0,0],'Localization',[1,0.5,0],'ZeroLine',[1,0,0],'EPhys',[0,0,0],'LickP1',[0,0,1],'LickP2',[0,0,1],'Cam64',[0,1,0]);
MV.VideoAvailable = any(strcmp(MV.Modules,MV.VideoName)) * P.ShowVideo;
MV.AudioAvailable = any(strcmp(MV.Modules,'AnalogIn')) * (P.ShowAudio | P.ShowSpec);
MV.LocalizationAvailable = any(strcmp(MV.Modules,'Localization'));
MV.EPhysAvailable = any(strcmp(MV.Modules,'EPhys')) * P.ShowEPhys;
MV.VocPosition = NaN;

MV.Trials = R.NIDAQ.Trials; MV.NTrials = length(MV.Trials); 
R.NIDAQ.SRAI = 1/diff(R.NIDAQ.Data(1).Data.Time(1:2));
if any(diff(MV.Trials)>1) error('Only consecutive trials are supported in MultiViewer at thie point.'); end

for iT =1:MV.NTrials
  L = double(size(R.NIDAQ.Data(iT).Data.Time,1));
  MV.DAQLengths(iT) = L/R.NIDAQ.SRAI;
end
% COLLECT START TIMES FOR DIFFERENT DATA SOURCES
MV.DAQStarts = [0,cumsum(MV.DAQLengths(1:end-1))];
if MV.EPhysAvailable
  for iS=1:length(R.NIDAQ.Data)
    MV.EPhysStarts(iS) = R.EPhys.Data(iS).Data.Time(1);
  end
end
% VIDEO
if P.ShowVideo && MV.VideoAvailable
  for iT=1:length(R.NIDAQ.Data)
    switch MV.VideoName
      case {'Thorcam','Scanimage','VideoCalcium'}
        MV.VideoStarts(iT) = R.NIDAQ.Data(iT).Data.Time(1); 
      otherwise
        if ~isempty(R.NIDAQ.Data(iT).CameraStartTime)
          MV.VideoStarts(iT) = R.NIDAQ.Data(iT).CameraStartTime(1);
        else
          MV.VideoStarts(iT) = NaN;
        end
    end
  end
end

% LOCALIZATION
if MV.LocalizationAvailable
  for iT=1:length(R.NIDAQ.Data)
    if ~isempty(R.NIDAQ.Data(iT).CameraStartTime)
      MV.LocalizationStarts(iT) = R.NIDAQ.Data(iT).LocalizationStartTimes(1);
    else
      MV.LocalizationStarts(iT) = NaN;
    end
  end
end

MV.TimeTotal = sum(MV.DAQLengths);
MV.StartTime = R.NIDAQ.Data(1).Data.Time(1,1);
if ~isempty(R.NIDAQ.Data(end).Data.Time)
  MV.StopTime = R.NIDAQ.Data(end).Data.Time(end,1);
else
  MV.StopTime = R.NIDAQ.Data(end-1).Data.Time(end,1);
end

if MV.P.ShowVideo
  if isfield(MV.SetupInfo.Video,'Mirror1') MV.P.Mirror1 = MV.SetupInfo.Video.Mirror1; end
  if isfield(MV.SetupInfo.Video,'Mirror2') MV.P.Mirror2 = MV.SetupInfo.Video.Mirror2; end
  if isfield(MV.SetupInfo.Video,'Transpose') MV.P.Transpose = MV.SetupInfo.Video.Transpose; end
  if isfield(MV.SetupInfo.Video,'CorrectDistortion') MV.P.CorrectDistortion = MV.SetupInfo.Video.CorrectDistortion; end
end

%% PREPARE FIGURE (COMPUTE SIZES)
SS = get(0,'ScreenSize');

figure(P.FIG); set(P.FIG,'Position',[50,400,SS(3:4)-[100,480]],'Toolbar','figure',...
  'Name',['MultiViewer ',R.General.Parameters.General.Identifier,' ',MV.SetupInfo.SetupName],'NumberTitle','off',...
  'DeleteFcn',{@MV_saveAnnotations});

DC = axesDivide(1,[1,0.1],[0.05,0.04,0.93,0.9],[ ],[0.15]);
if MV.VideoAvailable
  DCTop = axesDivide([0.5,0.5],1,DC{1},[0.05],[]);
else
  DCTop{2} =  DC{1};
end
YDivision = []; YSep = []; LastInd = 0; AudioInd = [];
if MV.AudioAvailable
  if MV.P.ShowAudio
    YDivision(end+1) = 0.2; 
    AudioInd = LastInd + 1;
    LastInd = AudioInd;
  end
  if MV.P.ShowSpec
    YDivision(end+1) = 0.2;
    YSep(end+1) = 0.1;
    AudioInd = [AudioInd,LastInd  + 1];
  end
  LastInd = AudioInd(end);
end
if MV.EPhysAvailable 
  YDivision(end+[1]) = [0.8];
  EPhysInd = LastInd+1; LastInd = EPhysInd; 
end
if MV.P.ShowSensors 
  YDivision(end+[1,2]) = [0.2];  
  SensorInd = LastInd + [1,2]; LastInd = SensorInd(end); 
end
DCTopRight = axesDivide(1,YDivision,DCTop{2},[],0.07);
DCBottom     = axesDivide([1.5,0.1,0.05,0.05,0.05,0.05],1,DC{2},[0.01],[]);

DS.Scrobble = DCBottom{1};
if MV.VideoAvailable  
  DS.Video                = DCTop{1}; 
  DS.Colorbar          = [DS.Video([1])+1.02*DS.Video([3]),DS.Video(2),0.01,DS.Video(4)];
end
if MV.AudioAvailable 
  DS.Spectrogram   = DCTopRight{AudioInd(1)}; 
  if MV.P.ShowAudio
    DS.Audio               = DCTopRight{(AudioInd(2))};
  end
end
if MV.P.ShowSensors
  DS.Triggers           = DCTopRight{SensorInd(1)};
  DS.Sensors            = DCTopRight{SensorInd(2)};
end
if MV.EPhysAvailable
  DS.EPhys                = DCTopRight{EPhysInd};
end
colormap(gray(256)); 

F = fieldnames(DS);
for i=1:length(F); MV.AH.(F{i}) = axes('Pos',DS.(F{i})); box on; end

%% CONTROLS
MV.CurrentTime = 0;  MV.CurrentSelInd = struct('left',1,'right',1,'center',1); MV.CurrentSelIndShift = 1; MV.Playing = 0;
MV.GUI.CurrentTime = uicontrol('style','edit','String',num2str(MV.CurrentTime),...
  'Units','normalized','Position',DCBottom{2},'Callback',{@MV_showData,'setedit'},'Tooltip','Current Time');
MV.GUI.StepSize = uicontrol('style','edit','String',num2str(MV.P.StepSize),...
  'Units','normalized','Position',DCBottom{3},'Tooltip','Step Size for Buttons');
MV.GUI.Help = uicontrol('style','pushbutton','String','Doc',...
  'Units','normalized','Position',DCBottom{4},'Tooltip','Show Documentation','Callback','doc MultiViewer');

ButtonStyles = {'pushbutton','pushbutton'};
switch MV.P.StepTo
  case 'Annotation'; StepCommand = 'StepAnno';
  case 'Vocalization'; StepCommand = 'StepVoc';
end
Commands = {'Step','Step';StepCommand,StepCommand};
Strings = {'<','>'};
Tooltips = {'Step Forward','Step Backward'};
W = [1,1,1]; Colors = {0.8*W,0.8*W};
for i=1:length(Strings)
  UH(i) = uicontrol('style',ButtonStyles{i},'String',Strings{i},'Units','normalized','Position',DCBottom{4+i},...
    'Callback',{@MV_showData,Commands{1,i},Strings{i}},'BackGroundColor',Colors{i},'Tooltip',Tooltips{i},'KeypressFcn',@LF_KeyPress);
  uicontrol('style',ButtonStyles{i},'String',Strings{i},'Units','normalized','Position',DCBottom{4+i} + [0,0.08,0,-0.02],...
    'Callback',{@MV_showData,Commands{2,i},Strings{i}},'BackGroundColor',0.8*Colors{i},'Tooltip',Tooltips{i},'KeypressFcn',@LF_KeyPress);
end

%% EXPORT CONTROLS
MV.GUI.StartTime = uicontrol('style','edit','String',num2str(MV.StartTime),...
  'Units','normalized','Position',DCBottom{2}./[1,1,1,1.5]+[0,1.1*DCBottom{2}(4),0,0],'Tooltip','Export Start Time');
MV.GUI.StopTime = uicontrol('style','edit','String',num2str(MV.StartTime),...
  'Units','normalized','Position',DCBottom{2}./[1,1,1,1.5]+[1.1*DCBottom{2}(3),1.1*DCBottom{2}(4),0,0],'Tooltip','Export Stop Time');

if MV.P.ShowSensors
  %% PREPARE DATA DISPLAY
  InputNames = {R.General.Paradigm.HW.Inputs.Name};
  MV.GUI.NIDAQ.SensorInd = find(~cellfun(@isempty,strfind(InputNames,'Pos')));
  axes(MV.AH.Sensors); set(MV.AH.Sensors,'ButtonDownFcn',{@MV_setTime}); hold on;
  MV.GUI.NIDAQ.SensorNames = InputNames(MV.GUI.NIDAQ.SensorInd);
  set(gca,'XLim',[-MV.P.Window/2,MV.P.Window/2],'YLim',[0,9]);
  if MV.P.ShowZeroLine;  plot([0,0],[-10,10],'Color',MV.Colors.ZeroLine); end
  MV.Colors.AniPosP1S1 = [0.5,0,0];
  MV.Colors.AniPosP1S2 = [1,0,0];
  MV.Colors.AniPosP1S3 = [1,0.3,0.3];
  MV.Colors.PlatPosP1 = [1,0.6,0.6];
  MV.Colors.AniPosP2S1 = [0,0,0.5];
  MV.Colors.AniPosP2S2 = [0,0,1];
  MV.Colors.AniPosP2S3 = [0.3,0.3,1];
  MV.Colors.PlatPosP2 = [0.6,0.6,1];
  MV.Colors.CamPosC1 = [0.3,1,0.3];
  
  for i = 1:length(MV.GUI.NIDAQ.SensorInd)
    cColor = MV.Colors.(MV.GUI.NIDAQ.SensorNames{i});
    MV.GUI.NIDAQ.SensorH(i) = plot(0,0,'Color',cColor,'Hittest','off');
    text(0.02+(i-1)*0.11,0.9,MV.GUI.NIDAQ.SensorNames{i},'Units','n','FontSize',8,'FontWeight','bold','Color',cColor,'Hittest','off');
  end
  ylabel('Voltage [V]');
  
  axes(MV.AH.Triggers); set(MV.AH.Triggers,'ButtonDownFcn',{@MV_setTime}); hold on;
  MV.GUI.NIDAQ.TriggerInd = setdiff([1:length(InputNames)],MV.GUI.NIDAQ.SensorInd);
  MV.GUI.NIDAQ.TriggerNames = InputNames(MV.GUI.NIDAQ.TriggerInd);
  MV.Colors.Trial = [1,0,0];
  MV.Colors.Recording = [1,0,0];
  MV.Colors.CamStart = [0,0,1];
  MV.Colors.CamTrigTo = [0,1,0];
  MV.Colors.CamTrigFrom = [0.5,1,0.5];
  MV.Colors.Camera = [0.5,1,0.5];
  MV.Colors.ThorcamTrigFrom = [0,0,1];
  MV.Colors.ThorcamTrigTo = [0,0.5,1];
  MV.Colors.ThorcamFrameValid = [0,0,0.5];
  MV.Colors.Feedback1 = [0,0,0];
  MV.Colors.Feedback2 = [0.5,0.5,0.5];
  MV.Colors.LED1 = [0.5,0.5,0.5];
  MV.Colors.Speaker = [0,1,1];  
  MV.Colors.Speaker = [0,0,0.5];
  MV.Colors.Reward = [0.5,1,0.5];
  MV.Colors.Opto = [0.5,1,0.5];
  MV.Colors.Lick = [0,0,1];
  MV.Colors.Vacuum = [0.5,1,0.5];
  MV.Colors.Frame = [0.5,1,0.5];

  NTriggers = length(MV.GUI.NIDAQ.TriggerInd);
  set(gca,'XLim',[-MV.P.Window/2,MV.P.Window/2],'YLim',[-0.1,2*NTriggers+0.1]);
  if MV.P.ShowZeroLine; plot([0,0],[-10,10],'Color',MV.Colors.ZeroLine); end
  for i = 1:NTriggers
    if isfield(MV.Colors,MV.GUI.NIDAQ.TriggerNames{i})
      cColor = MV.Colors.(MV.GUI.NIDAQ.TriggerNames{i});
    else
      cColor = [0,0,0];
    end
    MV.GUI.NIDAQ.TriggerH(i) = plot(0,0,'Color',cColor,'Hittest','off');
    text(i*0.1,0.9,[MV.GUI.NIDAQ.TriggerNames{i},' (',num2str(MV.GUI.NIDAQ.TriggerInd(i)),')'],'Units','n','FontSize',12,'FontWeight','bold','Color',cColor,'HitTest','off');
  end
  MV.GUI.NIDAQ.FrameH = plot(0,0,'.','Color',[0,0,1],'HitTest','off','MarkerSize',10);
  ylabel('Voltage [V]');
  xlabel('Time [s]');
  MV_showNIDAQ;
end

if MV.AudioAvailable
  %% PREPARE AUDIO DISPLAY
  R.AnalogIn.SRAI = R.General.Parameters.Setup.Audio.SRAI;
  if ~isfield(MV,'SetupInfo') MV.P.NAudioChannels = length(MV.SetupInfo.Audio.Channels); end
  if isempty(MV.P.NAudioChannels)
    R.AnalogIn.NChannels = size(R.AnalogIn.Data(1).Data.Analog,2);
  else
     R.AnalogIn.NChannels = MV.P.NAudioChannels;
  end
  SRAI = R.AnalogIn.SRAI;
  if isnan(MV.P.SDAudio)
    for i=1: R.AnalogIn.NChannels
      MV.GUI.Audio.AudioStd(i) = std(R.AnalogIn.Data(1).Data.Analog(1*SRAI:1.5*SRAI,i));
    end
  else
    MV.GUI.Audio.AudioStd = repmat(MV.P.SDAudio,1,R.AnalogIn.NChannels);
  end
  
  NWindow = 0.05*SRAI;
  [b,a] = butter(2,[20000/(SRAI/2)],'high');
  for i=1:length(R.AnalogIn.Data)
    if ~isempty(R.AnalogIn.Data(i).Data.Analog)
      R.AnalogIn.Data(i).Data.Analog = ...
        bsxfun(@times,...
        R.AnalogIn.Data(i).Data.Analog(:,1:R.AnalogIn.NChannels),...
        0.3./MV.GUI.Audio.AudioStd);
    end
    NSteps = length(R.AnalogIn.Data(i).Data.Analog);
    NSel = NWindow * floor(NSteps/NWindow);
    AnalogStd{i} = sum(abs(reshape(filter(b,a,R.AnalogIn.Data(i).Data.Analog(1:NSel)),NWindow,NSel/NWindow)))/NWindow;
    TimeStd{i} = R.AnalogIn.Data(i).Data.Time(NWindow/2:NWindow:NSel-NWindow/2);
  end
  MV.Colors.AnalogInTrace = {[1,0,0],[0,0,1],[0,1,0],[1,0.5,0]};
  R.AnalogIn.MicNames = MV.SetupInfo.Audio.MicNames;
  
  if MV.P.ShowAudio
    axes(MV.AH.Audio); set(MV.AH.Audio,'ButtonDownFcn',{@MV_setTime}); hold on;
    set(MV.AH.Audio,'YDir','normal','YLim',[-10,10],'XLim',[-MV.P.Window/2,MV.P.Window/2]);
    for i = 1:R.AnalogIn.NChannels
      MV.GUI.Audio.SoundH(i) = plot(MV.AH.Audio,0,0,'-','Color',MV.Colors.AnalogInTrace{i},'Hittest','off');
      text(0.01+(i-1)*0.02,0.07,R.AnalogIn.MicNames{i},'Units','n','FontSize',10,'FontWeight','bold','Color',MV.Colors.AnalogInTrace{i},'Hittest','off','Units','n');
    end
    if MV.P.ShowZeroLine; plot([0,0],[-10,10],'Color',MV.Colors.ZeroLine,'Hittest','off'); end
    ylabel('Voltage [V]');
end

  %% PREPARE SPECTROGRAM DISPLAY
  if MV.P.ShowSpec
    axes(MV.AH.Spectrogram); set(MV.AH.Spectrogram,'ButtonDownFcn',{@MV_setTime}); hold on;
    MV.GUI.Audio.SpectrogramH = imagesc(0,0,0,'Hittest','off');
    if MV.P.ShowZeroLine; plot([0,0],[0,R.AnalogIn.SRAI/2],'Color',MV.Colors.ZeroLine,'HitTest','off'); end
    
    R.AnalogIn.CorrTime = 0.0001;
    XTime = 1000*[-R.AnalogIn.CorrTime:1/SRAI:R.AnalogIn.CorrTime];
    MV.GUI.Audio.XCorrH = plot(XTime,10000*ones(size(XTime)),'g','Hittest','off');
    MV.GUI.Audio.DiffPosH = text(0.5,1.05,'','Color','r','FontSize',14,'FontWeight','bold','Horiz','center','HitTest','off','Units','n');
    set(MV.AH.Spectrogram,'YDir','normal','YLim',[0,125],'XLim',[-MV.P.Window/2,MV.P.Window/2]);
    ylabel('Freq. [kHz]');
    caxis([0,20]);
    MV.GUI.Audio.SpecChannel = 1;
    MV.GUI.Audio.Source = 'AnalogIn';
    MV.GUI.SpecChannel = text(0,0.2,R.AnalogIn.MicNames{MV.GUI.Audio.SpecChannel},...
     'Units','n','FontSize',14,'FontWeight','Bold','Color',MV.Colors.AnalogInTrace{MV.GUI.Audio.SpecChannel});
    
    %% ADD THE VOCALIZATIONS
    if MV.P.DetectVocs 
      MV_DetectVocs;
    else
      if isempty(MV.P.Vocs) % Check on Server
        cPar = R.General.Parameters.General;
        [~,Paths] = C_getDir('Animal',cPar.Animal,'Recording',cPar.Recording);
        File = [Paths.DB,'Results',filesep,'Vocalizations',filesep,'Vocalizations.mat'];
        if exist(File,'file')
          fprintf(['Loading previously collected Vocalizations from [ ',escapeMasker(File),' ]\n']);
          tmp = load(File); MV.P.Vocs = tmp.Vocs;
        end
      end
    end
    
    if ~isempty(MV.P.Vocs)
      Opts = {'.','HitTest','off','MarkerSize',14};
      MV.GUI.Audio.VocStarts = plot(zeros(size(MV.P.Vocs)),110*ones(size(MV.P.Vocs)),Opts{:},'Color','g');
      MV.GUI.Audio.VocStops = plot(zeros(size(MV.P.Vocs)),110*ones(size(MV.P.Vocs)),Opts{:},'Color','r');
    end
  end
  
  MV_showAudio;
end

%% PREPARE VIDEO DISPLAY
if MV.VideoAvailable
  axes(MV.AH.Video); hold on;
  R.(MV.VideoName).Dims = R.(MV.VideoName).Data(1).Data.Dims(1:2);   % Dimensions are typically 512x640
  MV.GUI.Video.NFrames = 0;
  for i=1:length(R.(MV.VideoName).Data)
    MV.GUI.Video.NFrames = MV.GUI.Video.NFrames + R.(MV.VideoName).Data(i).Data.Dims(end);
  end
  if isempty(MV.P.Annotations) MV.Annotations = cell(MV.GUI.Video.NFrames,1); 
    for iF=1:MV.P.RandomFrames 
      cFrame = ceil(MV.GUI.Video.NFrames*rand); if cFrame == 0 cFrame = 1; end
      MV.Annotations{cFrame} = NaN;
    end
  else MV.Annotations = MV.P.Annotations;
  end
  MV.CurrentFrame = 0;
  
  NBinSteps  = 11;
  for i=1:2
    Steps{i} = linspace(MV.SetupInfo.Video.ViewField(1,i),MV.SetupInfo.Video.ViewField(2,i),MV.SetupInfo.Video.Dims(i));
    Ticks{i} = HF_autobin(MV.SetupInfo.Video.ViewField(1,i),MV.SetupInfo.Video.ViewField(2,i),NBinSteps);
  end
  MV.GUI.Video.FrameH = imagesc(MV.SetupInfo.Video.XPos,MV.SetupInfo.Video.YPos,zeros(MV.SetupInfo.Video.Dims));
  set(MV.AH.Video,'YDir','normal');
  set(MV.GUI.Video.FrameH,'HitTest','off');  hold on;
 MV.GUI.Video.ClassH = text(0.02,0.02,'','Units','n','Color','r','FontW','b');
 MV.GUI.Video.ClassShiftH = text(0.98,0.02,'','Units','n','Color','r','FontW','b','Horiz','r');
 MV.GUI.Video.ClassC = text(0.42,0.02,'','Units','n','Color','r','FontW','b'); % added by G, display # annotated frames on screen
          
  set(MV.AH.Video,'XTick',Ticks{1},'YTick',Ticks{2},...
    'PlotBoxAspectRatio',[MV.SetupInfo.Video.Dims./max(MV.SetupInfo.Video.Dims),1],...
    'PlotBoxAspectRatioMode','manual',...
    'ButtonDownFcn',{@MV_mainFrame},'XLim',MV.SetupInfo.Video.ViewField(:,1),'YLim',MV.SetupInfo.Video.ViewField(:,2));
  
  if MV.P.ShowGrid
    for i=1:length(Ticks{1})
      plot3([Ticks{1}(i),Ticks{1}(i)],Ticks{2}([1,end]),[256,256],':','Color','w','HitTest','off');
    end
    for i=1:length(Ticks{2})
      plot3(Ticks{1}([1,end]),[Ticks{2}(i),Ticks{2}(i)],[256,256],':','Color','w','HitTest','off');
    end
   end
  MV.GUI.Video.VocLocation = plot3([-40],[-40],[256,256],'+','Color',[1,0.5,0],'Hittest','off');
  xlabel('X [m]'); ylabel('Y [m]');
  caxis([100,255]);
  MV.GUI.Video.TitleH = title(['Video']);
  MV.GUI.AllAnnotations = plot([-1000],[-1000],'.r'); 
  set(MV.GUI.AllAnnotations,'HitTest','off');
  
  % ADD THE MICROPHONES
  if MV.AudioAvailable
    Mics = MV.SetupInfo.Audio.MicNames;
    MicPos = MV.SetupInfo.Audio.MicrophonePositions;
    for iM=1:length(Mics)
      text(MicPos{iM}(1),MicPos{iM}(2),Mics{iM},'FontSize',8,'FontWeight','bold','Color',MV.Colors.AnalogInTrace{iM});
    end
  end
  
  % ADD THE LOCALIZATION
  if MV.LocalizationAvailable
    NX = length(R.Localization.X); NY = length(R.Localization.Y);  BaseFrame = ones(NY,NX);
    MV.GUI.Localization.FrameH = imagesc(R.Localization.X,R.Localization.Y,0*BaseFrame);
    MV.GUI.Localization.MaxH = plot(0,0,'y.','MarkerSize',20);
    %MV.GUI.Localization.DirH = plot(0,0,'r-');
    set(MV.GUI.Localization.FrameH,'AlphaData',0.5*BaseFrame);
    set(MV.GUI.Video.FrameH,'HitTest','off');
    R.Localization.Dims = R.Localization.Data(1).Data.Dims(1:2);  
    MV.Localization.Median = median(R.Localization.Data.Data.Frames(:));
    MV.Localization.Min = min(R.Localization.Data.Data.Frames(:));
    MV.Localization.Max = max(R.Localization.Data.Data.Frames(:));
    MV.Localization.Range = MV.Localization.Max - MV.Localization.Min;
  end
    
  %% COLORBAR & COLOR CONTROLS
  axes(MV.AH.Colorbar);
  MV.NColSteps = 256; colormap(gray(MV.NColSteps)); MV.Colormap =0;
  MV.CLimits = [MV.P.CMin,MV.P.CMax]; 
  MV.GUI.CIH = imagesc(1); set(MV.GUI.CIH,'Hittest','off');
  set(MV.AH.Colorbar,'ButtonDownFcn',{@MV_setCLimit,'setvisual'})
  MV_setCLimit;
  
  %% ADD ANNOTATIONS
  if P.LoadAnnotations
    switch class(R.General.Paradigm)
      case {'AudioVideoRecordingLogic'}
        MV_loadAnnotations([],[],'Manual');
    end
  end
end

%% PREPARE EPHYS DISPLAY
if MV.EPhysAvailable
  SRAI = R.EPhys.SRAI;
  NElectrodes = length(MV.P.Electrodes);
  Range = 200;
  axes(MV.AH.EPhys); set(MV.AH.EPhys,'ButtonDownFcn',{@MV_setTime}); hold on;
  set(MV.AH.EPhys,'YDir','normal','YLim',[0,Range*NElectrodes],'XLim',[-MV.P.Window/2,MV.P.Window/2],'YTick',[]);
  
  for iE=1:NElectrodes
    if mod(floor((iE-1)/16),2)==0 cColor = [0,0,0.5]; else cColor = [0.5,0,0]; end
    MV.GUI.EPhys.DataH(iE) = plot(MV.AH.EPhys,0,0,'-','Color',cColor,'Hittest','off');
    text(-1.01,Range*(iE-0.5),['El.',num2str(MV.P.Electrodes(iE))],...
      'FontSize',12,'FontWeight','bold','Color',cColor,'Hittest','off','Horiz','right');
  end
  plot([0,0],[0,(iE+1)*Range],'r');
  xlabel('Time [s]');
  
  if MV.P.ShowSpikes
    % TRY TO LOAD SPIKES
    fprintf('Loading Spikes...\n');
    MV.EPhys.Resp = C_getResponse('Animal',P.Animal,'Recording',P.Recording,'ResponseType','Singleunit','Electrodes',P.Electrodes,'Reference','all');
  
    iK=0;
    for iES=1:length(MV.EPhys.Resp) % LOOP OVER ELECTRODE SETS
      cResp = MV.EPhys.Resp(iES);
      if iscell(cResp)
          cResp = cResp{1};
      end
       for iPr=1:length(cResp) % LOOP OVER ELECTRODE SETS
           for iC=1:length(cResp(iPr).NSpikes) % LOOP OVER CELL
               cElectrode = cResp(iPr).BestElectrodeByCell(iC);
               cInd = find(cElectrode==MV.P.Electrodes);
               if ~isempty(cInd)
                   iK=iK+1; MV.GUI.EPhys.SpikesH(iK) = plot(0,0,'r.','MarkerSize',12);
               end
           end
       end
    end
  end
  
  MV.GUI.EPhys.Electrodes = MV.P.Electrodes;
  ylabel('Voltage [\muV]','Units','n','Position',[-0.05,0.5]);
end

%% CREATE SCROBBLER
% Plot Markers at the Trial ends (red bars)
axes(MV.AH.Scrobble); hold on;
plot3([MV.CurrentTime,MV.CurrentTime],[0,1],[1,1],'r','Hittest','off'); 
axis([MV.StartTime,MV.StopTime,0,1]);
% SHOW ANIMAL POSITION
if 0 % ~isempty(R.General.Paradigm.HW.AnimalSensorsPhys)
  [AnimalPositionVsTime , AnimalTime , AnimalPositions ] = MV_computeAnimagesc(AnimalTime,AnimalPositions,AnimalPositionVsTime,'HitTest','off');
  caxis([-1,1]);
end

for iT = 1:MV.NTrials
  if ~isempty(R.NIDAQ.Data(iT).TrialStartTime)
    MV.TrialStarts(iT) = double(R.NIDAQ.Data(iT).TrialStartTime);
    MV.TrialStops(iT) = double(R.NIDAQ.Data(iT).TrialStopTime);
    patch([MV.TrialStarts(iT),MV.TrialStops(iT),MV.TrialStops(iT),MV.TrialStarts(iT)],[0,0,1,1],[0.1,0.1,0.1,0.1],0,'FaceColor',[0.8,0.8,0.8],'HitTest','off','facealpha',0.5)
    text(MV.TrialStarts(iT),1.3,num2str(iT),'Color','k','Horiz','center','FontSize',8);
  end
end

% Plot Scrobbler bar
MV.GUI.CurrentLine = plot3(repmat(MV.CurrentTime,1,2), [0,1],[1,1],'Color',MV.Colors.ZeroLine,'HitTest','off','LineWidth',2);
if MV.VideoAvailable set(MV.AH.Video,'ZLim',[0,256]); end

% Plot horizontal bars in different colors to show when different modalities were recorded
for iT=1:MV.NTrials
  for iM=1:length(MV.Modules)
    cModule = MV.Modules{iM};
    cColor = MV.Colors.(cModule);
    ModPos = iM/(length(MV.Modules)+1);
    if iT == 1
      text(-1,ModPos,iM,[cModule,' '],'Horiz','right','Color',cColor,'FontWeight','bold');
    end
    if length(R.(cModule).Data)>=iT && ...
        ~isempty(R.(cModule).Data(iT).Data.Time)
        cTimes = R.(cModule).Data(iT).Data.Time([1,end],1);
          plot(cTimes,repmat(ModPos,1,2),'.-','Color',cColor,'LineWidth',1.5,'HitTest','off');
      switch cModule
        case MV.VideoName % CHECK FOR INTERNAL STOPS WHERE NO VIDEO WAS ACQUIRED (TO BE DONE)
          InterTimes = diff(R.(cModule).Data(iT).Data.Time(:,1));
          cDiffInds= [0;find(InterTimes>0.1);length(InterTimes)+1];
          for iP = 1:length(cDiffInds)-1
            cInds = [cDiffInds(iP)+1,cDiffInds(iP+1)];
            cTimes = R.(cModule).Data(iT).Data.Time(cInds,1);
            plot(cTimes,repmat(ModPos,1,2),'.-','Color',0.5*cColor,'LineWidth',1.5,'HitTest','off','MarkerSize',16);
          end
        case 'AnalogIn'
          if MV.AudioAvailable
            for i=1:length(TimeStd)
              cInd = intersect(find(AnalogStd{i}>0.25),find(AnalogStd{i}<1));
              plot(TimeStd{i}(cInd),repmat(ModPos,size(cInd)),'.','Color',cColor,'MarkerSize',16,'HitTest','off')
            end
          end
          if ~isempty(MV.P.Vocs)
            VocTimes = [MV.P.Vocs.Start];
            plot(VocTimes,repmat(ModPos,size(VocTimes)),'.','Color','g','MarkerSize',12,'HitTest','off')
          end
      end
    end
  end
end


% Activate Scrobbler to choose point in time
text(0,-0.25,'Time [s]  ','Units','n','FontWeight','bold','Horiz','right');
XTicks = C_autobin(MV.StartTime,MV.StopTime,10);
set(MV.AH.Scrobble,'ButtonDownFcn',{@MV_scrobblor},'XTick',XTicks,'YTick',[],'XGrid','on');
set(P.FIG,'WindowButtonUpFcn','global Scrobbling_ ; Scrobbling_ = 0;');
% ATTACH KEYPRESS FUNCTION LATE TO AVOID ACCIDENTAL CLICKING
set(P.FIG,'KeyPressFcn',{@LF_KeyPress});

if ~isempty(MV.P.Time); MV_setTime([],[],MV.P.Time); end

%% CALLBACK FUNCTIONS
function MV_mainFrame(O,E)
global MV R;

SelType = get(MV.FIG,'SelectionType');
NAnnotated = sum(cellfun(@isstruct,MV.Annotations)); 
set(MV.GUI.Video.ClassC,'String',NAnnotated); 

switch SelType
  case 'normal' % Mark left animal
    if MV.CurrentSelInd.left <= length(MV.P.SelTypes)
      MV_recordLocation(O,E,'left'); MV_setSelection( MV.CurrentSelInd.left+1,'left');
    end
  case 'alt' % Mark right animal
    if MV.CurrentSelInd.right <= length(MV.P.SelTypes)
      MV_recordLocation(O,E,'right'); MV_setSelection( MV.CurrentSelInd.right+1,'right');
    end
  case 'extend'
    if MV.CurrentSelInd.center <= length(MV.P.SelTypes)
      MV_recordLocation(O,E,'center'); MV_setSelection( MV.CurrentSelInd.center+1,'center');
    end
  case 'open' 
    if ~isempty(MV.Annotations{MV.CurrentFrame})
        MV.Annotations{MV.CurrentFrame}.Clicks = MV.Annotations{MV.CurrentFrame}.Clicks(1:end-1);
        MV_setSelection( max([MV.CurrentSelInd.left-1,0]),[],'left'); 
    end
    MV_zoomFrame(O,E); MV_showData([],[],'redraw')
end

function MV_setSelection(SelInd,ClickType,Modifier)
    global MV
    if ~SelInd cInd = 1; end
    if nargin<3
      if SelInd > length(MV.P.SelTypes) return; end
      %set(MV.GUI.Video.ClassH,'String',MV.P.SelTypes{SelInd}); 
      MV.CurrentSelInd.(ClickType) = SelInd;
    else
      switch Modifier
        case 'Shift'
          if SelInd > length(MV.P.SelTypesShift) SelInd = 1; end
          set(MV.GUI.Video.ClassShiftH,'String',MV.P.SelTypesShift{SelInd});
          MV.CurrentSelIndShift = SelInd;
      end
    end

    
function MV_zoomFrame(O,E)
global MV R;
Point = get(O,'CurrentPoint');
Point = Point(1,1:2);
XRange = diff(MV.SetupInfo.Video.XPos([1,end])); 
YRange = diff(MV.SetupInfo.Video.YPos([1,end])); 
set(MV.AH.Video,...
  'XLim',[Point(1) - XRange/4,Point(1)+XRange/4],...
  'YLim',[Point(2) - YRange/4,Point(2)+YRange/4]);

function MV_zoomOut(O,E)
global MV R
set(MV.AH.Video,'XLim',MV.SetupInfo.Video.ViewField(:,1)+[-0.01,0.01]','YLim',MV.SetupInfo.Video.ViewField(:,2)+[-0.01,0.01]');

function MV_recordLocation(O,E,ClickType)
global MV R;
if ~exist('ClickType','var') ClickType = ''; end
if MV.CurrentFrame
  switch ClickType
    case 'refresh'
    otherwise
      Point = get(O,'CurrentPoint');
      Point = Point(1,1:2);
      switch architecture
        case 'MAC'
          % LOCATION IS  CONSISTENTLY IMPRECISE, EMPIRICAL CORRECTION
          Point = Point + [-0.0002,0.0005];
      end
      [~,XPixel] = min(abs(MV.SetupInfo.Video.XPos - Point(1)));
      [~,YPixel] = min(abs(MV.SetupInfo.Video.YPos - Point(2)));
      fprintf(['X : ',num2str(Point(1)),'m  (',num2str(XPixel),' pixels)   Y : ',num2str(Point(2)),'m (',num2str(YPixel),' pixels)\n']);
      if isfield(MV,'CurrentFrame')
        if ~isstruct(MV.Annotations{MV.CurrentFrame}) & isnan(MV.Annotations{MV.CurrentFrame}) 
          MV.Annotations{MV.CurrentFrame} = []; end
        if isempty( MV.Annotations{MV.CurrentFrame}) cClickInd = 1; 
          MV.CurrentSelInd = struct('left',1,'right',1,'center',1);
        else cClickInd = length( MV.Annotations{MV.CurrentFrame}.Clicks)+1; end
        cData = struct('X',Point(1),'Y',Point(2),'XPixel',XPixel,'YPixel',YPixel,...
          'Time',MV.CurrentTime,'Frame',MV.CurrentFrame,...
          'Type',MV.P.SelTypes{MV.CurrentSelInd.(ClickType)},...
          'TypeShift',MV.P.SelTypesShift{MV.CurrentSelIndShift},...
          'VocLocation',MV.VocPosition,'ClickType',ClickType); 
        if isempty(MV.Annotations{MV.CurrentFrame}) || isempty(MV.Annotations{MV.CurrentFrame}.Clicks)
          MV.Annotations{MV.CurrentFrame}.Clicks = cData;
        else
          MV.Annotations{MV.CurrentFrame}.Clicks(cClickInd) = cData;
        end
        
      end
  end
  XData = find(~cellfun(@isempty,MV.Annotations));
  set(MV.GUI.AllAnnotations,'XData',XData,'YData',zeros(size(XData)) + 1);
  MV_showVideo;

end

%% CALLBACK FUNCTIONS
function LF_KeyPress(handle,event,FID)

global MV R
CC = get(MV.FIG,'CurrentCharacter');
if ~isempty(CC)
  switch int8(CC)
    case 28; MV_showData([],[],'step','<');  % left
    case 29; MV_showData([],[],'step','>'); % right
    case 32;  % space = play
      MV_showData([],[],'play','>');
    case 100; %d Delete Annotations for current Frame
      MV.Annotations{MV.CurrentFrame} = [];
      MV_recordLocation([],[],'refresh');
      MV.CurrentSelInd = struct('left',1,'right',1,'center',1);
      MV_showVideo;
    case 102; % f : jump to frame
      MV_gotoFrame;
    case 117; %u : undo last Annotation
      if ~isempty(MV.Annotations{MV.CurrentFrame}) & ~isempty(MV.Annotations{MV.CurrentFrame}.Clicks)
        LastClick = MV.Annotations{MV.CurrentFrame}.Clicks(end);
        MV.Annotations{MV.CurrentFrame}.Clicks = MV.Annotations{MV.CurrentFrame}.Clicks(1:end-1);
        if isempty(MV.Annotations{MV.CurrentFrame}.Clicks) MV.Annotations{MV.CurrentFrame} = []; end
        MV.CurrentSelInd.(LastClick.ClickType) = max([MV.CurrentSelInd.(LastClick.ClickType)-1,1]);
        MV_recordLocation([],[],'refresh');
        MV_showVideo;
      end
    case 113; % q : jump to previous annotation / 
      MV_showData([],[],'stepvoc','<');
    case 119; % w : jump to next annotation / 
      MV_showData([],[],'stepvoc','>');
    case 101; % e : export different kinds of data
      MV_exportData; 
    case 83; % S : change Audio channels
      switch MV.GUI.Audio.Source
        case 'AnalogIn'
          if MV.GUI.Audio.SpecChannel<R.AnalogIn.NChannels
            MV.GUI.Audio.SpecChannel = MV.GUI.Audio.SpecChannel + 1;
          else
            MV.GUI.Audio.Source = 'Cam64';
            MV.GUI.Audio.SpecChannel = NaN;
          end
        case 'Cam64'
          MV.GUI.Audio.Source = 'AnalogIn'; MV.GUI.Audio.SpecChannel=1;
      end
      if MV.GUI.Audio.SpecChannel==0  MV.GUI.Audio.SpecChannel = R.AnalogIn.NChannels; end
      MV_showData([],[],'redraw',''); 
    case 115; % s save annotations before quitting  
      MV_saveAnnotations; 
    case {43,61}; % =/+
      CLim = get(MV.AH.Spectrogram,'CLim'); set(MV.AH.Spectrogram,'CLim',1.2*CLim);
    case {45,95}; % -/_
      CLim = get(MV.AH.Spectrogram,'CLim'); set(MV.AH.Spectrogram,'CLim',0.8*CLim);
    case {48,49,50,51,52,53,54,55,56,57}; % 1,2,3,4,5 ... Selections, e.g. for Animal Tracking
      cInd = str2num(CC); MV_setSelection(cInd);
    case {33,64,35,36,37,94,38,42,40,41}; % Shift+1,2,3,4,5 ... Selections, e.g. for USV classification
      cInd = find(CC==[33,64,35,36,37,94,38,42,40,41]); MV_setSelection(cInd,'Shift');        
    case 111; % o open annotation
      MV_loadAnnotations;
    case 110; %n count annotated frames
      NAnnotated = sum(cellfun(@isstruct,MV.Annotations));
      fprintf('Number of Annotated Frames : %d\n',NAnnotated);
    case 112; % p show localizations
      MV_showAVLocalization;
    case 122; %z zoom out again
      MV_zoomOut;
  end
end
  
function MV_showData(O,E,Command,iFrame)
global MV R

switch lower(Command)
  case 'set'; MV.CurrentTime = iFrame;
  case 'setedit'; MV.CurrentTime= str2num(get(O,'String'));
  case 'setvisual'; MV.CurrentTime = get(O,'CurrentPoint');MV.CurrentTime = round(MV.CurrentTime(1,1));
  case 'step';   
    MV.P.StepSize = str2num(get(MV.GUI.StepSize,'String'));
    switch iFrame
      case '>'; MV.CurrentTime = MV.CurrentTime + MV.P.StepSize; 
      case '<'; MV.CurrentTime = MV.CurrentTime - MV.P.StepSize; 
    end
  case 'stepanno'; % Step to the next annotation
    if MV.VideoAvailable
      AnnotatedFrames = find(~cellfun(@isempty,MV.Annotations));
      if isempty(AnnotatedFrames) disp('No Annotations found!'); return; end
      if ~isempty(MV.CurrentFrame)
        FrameDist = AnnotatedFrames - MV.CurrentFrame;
        switch iFrame
          case '>';  FrameDist(FrameDist <= 0) = inf; [~,Pos] = min(FrameDist);
          case '<';  FrameDist(FrameDist >= 0) = -inf; [~,Pos] = max(FrameDist);
        end
        NewFrame = AnnotatedFrames(Pos);
      else
        NewFrame = AnnotatedFrames(1);
      end
      MV.CurrentTime = R.(MV.VideoName).Data.Data.Time(NewFrame);
    end
  case 'stepvoc'
    if MV.AudioAvailable
      if isempty(MV.P.Vocs) disp('No Vocalizations Loaded!'); return; end
      VocTimes = ([MV.P.Vocs.Start] + [MV.P.Vocs.Stop])/2;
      DiffTimes = VocTimes - MV.CurrentTime; DiffTime = [];
      switch iFrame
        case '>';  DiffTimes = DiffTimes(DiffTimes>0); if ~isempty(DiffTimes) DiffTime = DiffTimes(1); end
        case '<';  DiffTimes = DiffTimes(DiffTimes<0); if ~isempty(DiffTimes) DiffTime = DiffTimes(end); end
      end
      if isempty(DiffTime) disp('No Vocalization found in this direction!'); return; end
      MV.CurrentTime = MV.CurrentTime + DiffTime;
    end
  case 'redraw';
  case 'play';
    MV.Playing = ~MV.Playing;
    while MV.CurrentTime<MV.StopTime & MV.Playing
      MV_showData(O,E,'step','>');
      drawnow;
    end
end

MV.CurrentTime = min([MV.StopTime,MV.CurrentTime]);
MV.CurrentTime= max([MV.StartTime,MV.CurrentTime]);

set(MV.GUI.CurrentTime,'String',num2str(MV.CurrentTime));
if isfield(MV.GUI,'CurrentLine')
  set(MV.GUI.CurrentLine,'XData',[MV.CurrentTime,MV.CurrentTime]);
end

if MV.P.ShowSensors MV_showNIDAQ; end % Shows NIDAQ ( SENSORS and TRIGGERS )
if MV.AudioAvailable MV_showAudio; end % Shows AnalogIn (SoundPressure and Spectrogram) 
if MV.VideoAvailable MV_showVideo; end % Shows Video Data (FRAME)
if MV.LocalizationAvailable; MV_showLocalization; end % Shows LocalizationData (FRAME from Analysis of Cam64)
if MV.EPhysAvailable MV_showEPhys; end % Shows EPhys Data (Analog Data)

function MV_exportData;
  global MV R;
  
  % SHOW OPTIONS
  DataType = input('\nSelect Type of Data to Export: \n [A] Audio \n [V] Video \n [G] GUI\n ','s');
  TimeSelection = input('\nSelect Time Range To Export: \n [C] Current \n [S] Start to Stop Range\n [T] Tracked Frames\n','s');
  switch lower(TimeSelection)
    case 'c'; TimeSelection = 'Current';
    case 's'; TimeSelection = 'Range';
    case 't'; TimeSelection  = 'Tracked';
  end
  
  % GET FILENAME
  Filename = input('\nSpecify Filename [without extension, sent to Desktop]: ','s');
  switch architecture
    case 'PCWIN';
      [~,S]  =system('echo %HOMEPATH%');
      Dir = ['C:',S(1:end-1),filesep,'Desktop',filesep];
    case 'MAC';
      Dir = ['~/Desktop/'];
  end
  Filename = [Dir,Filename];
  
  % CALL SPECIFIC FUNCTIONS
  switch upper(DataType)
    case 'A'; MV_exportAudio(TimeSelection,Filename);
    case 'V';  MV_exportVideo(TimeSelection,Filename);
    case 'G'; MV_exportGUI(TimeSelection,Filename);
    otherwise fprintf('DataType not implemented.');
  end
  
function MV_exportAudio(TimeSelection,Filename)
  global MV R;
  if ~MV.AudioAvailable; fprintf('No Audio available.'); return; end
  switch TimeSelection
    case 'Current';
      fprintf('Exporting a single sound pressure point does not make sense.');
      return;
    case 'Range';
      StartTime =  str2num(get(MV.GUI.StartTime,'String'));
      StopTime =  str2num(get(MV.GUI.StopTime,'String'));
      cTrial = find(MV.DAQStarts<StartTime,1,'last');
      cTime = R.AnalogIn.Data(cTrial).Data.Time;
      cInd = find(cTime>=StartTime & cTime<=StopTime);
      cData = R.AnalogIn.Data(cTrial).Data.Analog(cInd,:);
    case 'Tracked';
      warning('Not implemented yet, ignoring.');
  end
  Filename = [Filename,'.mat'];
  fprintf(['Saving to ',escapeMasker(Filename),'\n']);
  save(Filename,'cData');
  
function MV_exportVideo(TimeSelection,Filename)
  global MV R;
  if ~MV.VideoAvailable; fprintf('No Videoavailable.'); return; end
  Format = input('\nChoose Format: [M]atlab, M[P]EG-4, [A]VI or [J]PG :','s');
  switch TimeSelection
    case 'Current';
      cData = get(MV.GUI.Video.FrameH,'CData');
    case 'Range';
      StartTime =  str2num(get(MV.GUI.StartTime,'String'));
      StopTime =  str2num(get(MV.GUI.StopTime,'String'));
      cTrial = find(MV.DAQStarts<StartTime,1,'last');
      cTime = R.(MV.VideoName).Data(cTrial).Data.Time;
      Frames = find(cTime>=StartTime & cTime<=StopTime);
      cData = MV_getVideoData(1,Frames);
    case 'Tracked';
      Frames = MV_getTrackedFrames;
      cData = MV_getVideoData(1,Frames);
  end
  fprintf(['Saving to ',escapeMasker(Filename),'\n']);
  
  switch upper(Format)
    case 'J'; 
      try delete([Filename,filesep,'*.jpg']); rmdir([Filename,filesep]); end
      mkdirAll([Filename,filesep]);
      for iF=1:size(cData,4)
        imwrite(repmat(squeeze(cData(:,:,1,iF)),[1,1,3,1]),...
          [Filename,filesep,'img',sprintf('%06d',Frames(iF)),'.jpg'],'jpg',...
          'Comment',['Frame ',num2str(Frames(iF))]);
      end
    case 'M'; save(Filename,'cData');
    case {'A','P'}; 
      switch upper(Format)
        case 'A'; MovieFormat = 'Motion JPEG AVI';
        case 'P'; MovieFormat = 'MPEG-4';
      end
      V = VideoWriter(Filename,MovieFormat);
      V.Quality = 50; open(V); writeVideo(V,permute(cData,[2,1,3,4])); close(V);
    otherwise fprintf('Format not known, assuming Matlab\n'); save(Filename,'cData');
  end
  
function MV_exportGUI(TimeSelection,Filename)

  global MV R;
  Format = input('\nChoose Format: M[P]EG-4 or [A]VI :','s');
  switch upper(Format)
    case 'A'; MovieFormat = 'Motion JPEG AVI';
    case 'P'; MovieFormat = 'MPEG-4';
  end
  
  switch TimeSelection
    case 'Current'
    case 'Range'
      StartTime =  str2num(get(MV.GUI.StartTime,'String'));
      StopTime =  str2num(get(MV.GUI.StopTime,'String'));
      StepSize =  str2num(get(MV.GUI.StepSize,'String'));
      
      %Dir = [getenv('USERPROFILE'),filesep,'Desktop'];
      V = VideoWriter(Filename,MovieFormat);
      V.open; CR = [];
      FigureSize = get(MV.FIG,'Position');
      Rectangle = [0,100,FigureSize(3),FigureSize(4)-100];
      for cTime=StartTime:StepSize:StopTime
        Range= fprintf([CR,'Time ',num2str(cTime)]);
        Range = Range-length(CR)/2;
        MV_showData([],[],'set',cTime);
        F = getframe(MV.FIG,Rectangle);
        V.writeVideo(F);
        CR = repmat('\b',1,Range);
      end
      V.close; fprintf('\n');  
  end
 
function [cTime,cData] = MV_getCurrentData(Module)
  global MV R
  StartTime = MV.CurrentTime-MV.P.Window/2;
  StopTime = MV.CurrentTime+MV.P.Window/2;
  switch Module
    case {'NIDAQ','AnalogIn'}; Starts = MV.DAQStarts; Field = 'Analog'; GetData = 1; 
    case MV.VideoName; Starts = MV.DAQStarts; Field = 'Frames'; GetData = 0;
    case 'EPhys'; Starts = MV.EPhysStarts; Field = 'Spike'; GetData = 1;
    case 'Cam64'; Starts = MV.DAQStarts; Field = 'Analog'; GetData = 1;
  end
 
  if ~GetData cData = []; end
  
  TrialStart = find(Starts<StartTime,1,'last');
  TrialStop = find(StopTime>Starts,1,'last');
  FillZerosStart = isempty(TrialStart);
  FillZerosStop = isempty(TrialStop);
   
  if isempty(TrialStart) && isempty(TrialStop)
    cTime = []; cData = []; return;
  end
   % FIND THE STARTING AND STOPPING POINT FOR THE PLOTTING (USES COARSE
   % SEARCH FIRST FOR SPEED UP);
   SearchDurationStep = 0.1; % Seconds
   cSR = 1/diff(R.(Module).Data(1).Data.Time(1:2,1));
   SearchStep = round(SearchDurationStep * cSR);
   if FillZerosStart;     cIndStart = [];
   else
     cIndStartCoarse = find(R.(Module).Data(TrialStart).Data.Time(SearchStep:SearchStep:end,1) >= StartTime,1,'first');
     if ~isempty(cIndStartCoarse)
       cInds = cIndStartCoarse*SearchStep+[-SearchStep+1:0];
       cIndStartFine = find(R.(Module).Data(TrialStart).Data.Time(cInds,1) >= StartTime,1,'first');
       cIndStart = (cIndStartCoarse-1)*SearchStep + cIndStartFine;
     else 
       cIndStart = [];
     end
   end
   
   if FillZerosStop;      cIndStop = [];
   else
     cIndStopCoarse = find(R.(Module).Data(TrialStop).Data.Time(SearchStep:SearchStep:end,1) <= StopTime,1,'last');
     if ~isempty(cIndStopCoarse)
       cInds = [(cIndStopCoarse-1)*SearchStep+1  : min((cIndStopCoarse+1)*SearchStep,length(R.(Module).Data(TrialStop).Data.Time))];
       cIndStopFine = find(R.(Module).Data(TrialStop).Data.Time(cInds,1) <= StopTime,1,'last');
       cIndStop = (cIndStopCoarse-1)*SearchStep + cIndStopFine;
     else 
       cIndStop = [];
     end
    end

   if TrialStart == TrialStop
     Inds = [cIndStart:cIndStop];
     cTime = R.(Module).Data(TrialStart).Data.Time(Inds,1);
     if GetData 
       switch Module
         case 'Cam64'; 
           switch class(R.(Module).Data(TrialStart).Data.(Field){1})
             case 'memmapfile'; cData = R.(Module).Data(TrialStart).Data.(Field){1}.Data.D(:,Inds)';
             otherwise  cData = R.(Module).Data(TrialStart).Data.(Field){1}(Inds,:); 
           end
           
          otherwise cData = R.(Module).Data(TrialStart).Data.(Field)(Inds,:); 
       end
     end
   else
     if FillZerosStart
       cTime = R.(Module).Data(TrialStop).Data.Time(1:cIndStop,1);
       if GetData cData = R.(Module).Data(TrialStop).Data.(Field)(1:cIndStop,:); end
     elseif FillZerosStop
       cTime = R.(Module).Data(TrialStart).Data.Time(cIndStart:end,1);
       if GetData cData = R.(Module).Data(TrialStart).Data.(Field)(cIndStart:end,:); end
     else
       cTime = ...
         [R.(Module).Data(TrialStart).Data.Time(cIndStart:end,1);...
         R.(Module).Data(TrialStop).Data.Time(1:cIndStop,1)];
       if GetData
         cData = ...
           [R.(Module).Data(TrialStart).Data.(Field)(cIndStart:end,:);...
           R.(Module).Data(TrialStop).Data.(Field)(1:cIndStop,:)];
       end
     end
   end
   % find_halfspace_mex finds the next index below the closest index (from below)
   % cIndStart = find(R.(Module).Times > MV.CurrentTime-MV.P.Window/2,1,'first');
   %   IndStart = find_halfspace_mex(R.(Module).Times,MV.CurrentTime-MV.P.Window/2);
   % cIndStop = find(R.(Module).Times < MV.CurrentTime+MV.P.Window/2,1,'last');
   %   IndStop = find_halfspace_mex(MV.(Module).Times,MV.CurrentTime+MV.P.Window/2);
   if GetData && size(cTime,1) ~= size(cData,1) keyboard; end   
  
function cData = MV_getVideoData(cTrial,Frames)
  global MV R
  switch class(R.(MV.VideoName).Data(cTrial).Data.Frames)
    case 'memmapfile'; cData = R.(MV.VideoName).Data(cTrial).Data.Frames.Data.D(:,:,:,Frames);
    otherwise cData = R.(MV.VideoName).Data(cTrial).Data.Frames(:,:,:,Frames);
  end
     
function MV_showNIDAQ
   global MV R;
 
   [cTime,cData] = MV_getCurrentData('NIDAQ');
   YStep = 2;
   % FILL THE SENSOR PLOT
   for i=1:length(MV.GUI.NIDAQ.SensorInd)
     set(MV.GUI.NIDAQ.SensorH(i),'YData',cData(:,MV.GUI.NIDAQ.SensorInd(i)) +3*floor((i-1)/4),'XData',cTime-MV.CurrentTime);
   end

   % FILL THE TRIGGER PLOT
   for i=1:length(MV.GUI.NIDAQ.TriggerInd)
     set(MV.GUI.NIDAQ.TriggerH(i),'YData',cData(1:2:end,MV.GUI.NIDAQ.TriggerInd(i))/5+YStep*(i-1),'XData',cTime(1:2:end)-MV.CurrentTime);
   end
   
   % SHOW THE CHOSEN FRAME TIMES AS DOTS ABOVE 
   if MV.VideoAvailable
     cTime = MV_getCurrentData(MV.VideoName);
     switch MV.VideoName
       case 'Thorcam'; TrigName = 'ThorcamTrigFrom';
       case 'Scanimage'; TrigName = 'ScanimageTrigFrom';
       case 'VideoCalcium'; TrigName = 'VideoCalciumTrigTo';  
       otherwise TrigName = 'CamTrigFrom';
     end
     CamTriggers = YStep*(find( strcmp(MV.GUI.NIDAQ.TriggerNames,TrigName) )-1);
     set(MV.GUI.NIDAQ.FrameH,'YData',CamTriggers*ones(size(cTime)),'XData',cTime - MV.CurrentTime);
   end
   
function MV_showEPhys
  global MV R
  
  [cTime,cData] = MV_getCurrentData('EPhys');
  
  if ~isempty(cTime)
    % FILL THE ELECTRODE PLOT
    MV.GUI.EPhys.YCenters = 200*[0.5:length(MV.P.Electrodes)-0.5];
    for iE=1:length(MV.P.Electrodes)
      set(MV.GUI.EPhys.DataH(iE),'YData',cData(:,iE)+MV.GUI.EPhys.YCenters(iE),'XData',cTime-MV.CurrentTime);
    end
    
    % SHOW THE SPIKE TIMES
    if MV.P.ShowSpikes
      iK=0;
      for iES=1:length(MV.EPhys.Resp) % LOOP OVER ELECTRODE SETS
        cResp = MV.EPhys.Resp{1}(iES);
        if iscell(cResp); cResp = cResp{1}; end
        for iC=1:length(cResp.NSpikes) % LOOP OVER CELL
          cElectrode = cResp.BestElectrodeByCell(iC);
          cElInd = find(cElectrode == MV.P.Electrodes,1);
          if ~isempty(cElInd)
            iK = iK+1; STs = cResp.Time{iC};
            cInd = find((STs > MV.CurrentTime-1) & (STs < MV.CurrentTime+1));
            set(MV.GUI.EPhys.SpikesH(iK),'YData', MV.GUI.EPhys.YCenters(cElInd)*ones(size(cInd)),'XData',STs(cInd)-MV.CurrentTime);
          end
        end
      end
    end
  end
  
function MV_showAudio
  global MV R;
  SRAI = R.General.Parameters.Setup.Audio.SRAI;
  NFFT = 512;
  
  
  [cTime,cData] = MV_getCurrentData('AnalogIn');
  assert(size(cTime,1) == size(cData,1),'Lengths of Data and Time Vector do not Match for AnalogIn.');
  
  if ~isempty(cTime)
    % FILL THE TRIGGER PLOT
    if MV.P.ShowAudio
      if R.AnalogIn.NChannels > 1 ShiftRange = 3; else ShiftRange = 1; end
      for i=1:R.AnalogIn.NChannels
        Shift = i*ShiftRange/R.AnalogIn.NChannels - ShiftRange/2;
        set(MV.GUI.Audio.SoundH(i),'YData',cData(1:end,i) + Shift,'XData',cTime-MV.CurrentTime);
      end
    end
    
    % ALLOW FOR SHOWING SORAMA DATA
    switch MV.GUI.Audio.Source
      case 'Cam64';    [cTimeS,cDataS] = MV_getCurrentData('Cam64'); 
        cDataS = double(cDataS); cDataS = cDataS*(1/std(cDataS(:)));
      case 'AnalogIn'; cTimeS = cTime; cDataS = cData(:,MV.GUI.Audio.SpecChannel);
    end
    
    NChannels = size(cDataS,2);
    for iC=1:NChannels 
      [S(:,:,iC),F,T] = HF_specgram(cDataS(:,iC),NFFT,SRAI,[0,125000],NFFT/2,0,1);
      if iC==1 & NChannels>1; S(end,end,NChannels) = 0; end
    end
    if isequal(MV.GUI.Audio.Source,'Cam64'); S = mean(S,3); end

    set(MV.GUI.Audio.SpectrogramH,'CData',abs(S),'XData',T - (MV.CurrentTime-cTime(1)),'YData',F/1000);
    
  else
    if MV.P.ShowAudio
      for i=1:R.AnalogIn.NChannels
        set(MV.GUI.Audio.SoundH(i),'YData',[],'XData',[]);
      end
    end
    set(MV.GUI.Audio.SpectrogramH,'CData',NaN);
  end
  
  % SET LABEL
  switch MV.GUI.Audio.Source
    case 'AnalogIn'
      cString = R.AnalogIn.MicNames{MV.GUI.Audio.SpecChannel}; 
      cColor = (MV.Colors.AnalogInTrace{MV.GUI.Audio.SpecChannel}+1)/2;
    case 'Cam64'; cString = 'Cam64'; cColor = [1,1,1];
  end
  set(MV.GUI.SpecChannel,'String',cString,'Color',cColor);

  % PLOT VOCALIZATION
  if ~isempty(MV.P.Vocs)
    RelStarts = [MV.P.Vocs.Start] - MV.CurrentTime;
    RelStops = [MV.P.Vocs.Stop] - MV.CurrentTime;
    set(MV.GUI.Audio.VocStarts,'XData',RelStarts);
    set(MV.GUI.Audio.VocStops,'XData',RelStops);
    % PLOT PRECOMPUTED LOCATION
    RelTimes = ([MV.P.Vocs.Start] + [MV.P.Vocs.Stop])/2 - MV.CurrentTime;
    [MIN,cInd] = min(abs(RelTimes));
    if isfield(MV.P.Vocs,'Location')
      cLocation = MV.P.Vocs(cInd).Location;
      if isfield(MV.GUI,'Video')
        if MIN<0.1 && ~isempty(cLocation)
          set(MV.GUI.Video.VocLocation,'XData',cLocation(1),'YData',cLocation(2),'Visible',1);
        else
          set(MV.GUI.Video.VocLocation,'XData',NaN,'YData',NaN,'Visible',0);
        end
      end
    end
  end
  
  % PLOT USV LOCATION BASED ON AUDIO
  if MV.P.ComputeLocation
    if ~isempty(cTime)
      ccTime = cTime - MV.CurrentTime;
      Ind = abs(ccTime) < 0.03;
      ccData = cData(Ind,:);
      for iC=1:R.AnalogIn.NChannels Sounds{iC} = ccData(:,iC); end
      if size(ccData,1) > 10
        LocMethod = 'Spectral';
        RLoc = VocLocalizer('Sounds',Sounds,'HighPass',35000,'SR',SRAI,'CorrMethod',{'EWGCC'},...
          'MicrophonePositions',MV.SetupInfo.Audio.MicrophonePositions,'SourceHeight',MV.SetupInfo.Video.Offset{1}(3),...
          'CenterShift',MV.SetupInfo.Audio.CenterShift,...
          'EstimationMethod',[],'MaxMethod','Max','IntersectionMethod','Map',...
          'FIG',1000);
        MV.VocPosition = RLoc.StimPosMean;
        cSCorr = 125*RLoc.SCorrAll{1,end}/(3*max(RLoc.SCorrAll{1,end}));
        set(MV.GUI.Audio.XCorrH,'YData',cSCorr,'XData',1000*RLoc.CorrTime);
        CertaintyColor = [min(100*mean(RLoc.StimPosSD)/2,1),0,0];
        set(MV.GUI.Audio.DiffPosH,'String',sprintf('%0.3fm  ',RLoc.StimPosMean(:)),'Color',CertaintyColor);
        set(MV.GUI.Video.VocLocation,'XData',RLoc.StimPosMean(1),'YData',RLoc.StimPosMean(2),'Visible',MV.VideoAvailable);
      end
    end
    MV.P.ComputeLocation = 0;
  end
  
function MV_showVideo
  global MV R;
  
  Trial = find(MV.CurrentTime>=MV.VideoStarts,1,'last');
  if ~isempty(Trial) && Trial > 0
    [~,cFrame] = min(abs(R.(MV.VideoName).Data(Trial).Data.Time(:,1) - MV.CurrentTime));
    % cFrame = find(R.(MV.VideoName).Data(Trial).Data.Time(:,1) >= MV.CurrentTime,1,'first'); % Previous code for finding the frame, which causes a slight shift of the previous tracking results
    cTime = R.(MV.VideoName).Data(Trial).Data.Time(cFrame,1);
    MV.CurrentFrame = cFrame;
    if ~isempty(cFrame) && abs(cTime - MV.CurrentTime) <0.1 % if within 100 ms of an existing frame
      switch class(R.(MV.VideoName).Data(Trial).Data.Frames)
        case 'memmapfile'
          cData = squeeze(R.(MV.VideoName).Data(Trial).Data.Frames.Data.D(:,:,1,cFrame))'; % Transpose to render image  in XY as dims 1 & 2
        otherwise
          cData = squeeze(R.(MV.VideoName).Data(Trial).Data.Frames(:,:,1,cFrame))'; % Transpose to render image  in XY as dims 1 & 2
      end
    else
      cData = zeros(R.(MV.VideoName).Dims);
      cFrame = NaN;
    end
    
    if isfield(MV.P,'CorrectDistortion') & MV.P.CorrectDistortion > 0
      % Original Problem: Image distorted by the camera lens and position/angle of camera
      % Additional Problem: 
      % - DLCTracking performed on uncorrected video images
      % - Manual Tracking performed on corrected images in MultiViewer using Strategy 1 (below)
      % Solution: Translate Markers/Videos into each other by correction algorithms.
      % Terminology: 
      % - Distorted : Uncorrect image or marker
      % - Undistorted : Corrected image or marker
      % 
      % Strategy 1:  was frame-by-frame using the function C_correctImageDistortion (original, written by B)
      % This approach uses a default parameter of 0.4
      %
      % For DLC or Manually tracked Markers on distorted images images:
      % - C_correctImageDistortion_V1_1(Data,'Mode','Markers') % where Data is a struct
      % For Manually tracked Markers on images that were undistorted with C_correctImageDistortion in MV:
      % - C_correctImageDistortion_V1_1(Data,'Mode','Markers') % where Data is a Matrix (3D, Bodypart X Dimension X Frame)
      % - see correctMVAnnotations for how it is called (including inversion of Y and change of Data format) 
      % 
      % For future recordings, Videos should be corrected before tracking in DLC/Manual:
      % The call for a per-video correction:
      % - C_correctImageDistortion_V1_1('Mode',Image) 
      % where the parameters need to be adapted to the specific set of recordings
     %  - For the original Cam64 recordings, the parameters are:
     %    'RadialStrength',.33,'PerspectiveStrength',[1.5,5.2],'ZoomX',1,'ZoomY',1.009,'TranslationX',1.5,'Platform',[127,619,73,440],'Axes',[400,370]
     % - For future or other old recordings, other parameter will have to be used.
      % VideoDistortionCorrection shows how to call this for a set of videos
      
      %cData = C_correctImageDistortion_V1_1(cData,'Mode','Image');
      cData = C_correctImageDistortion(cData,MV.P.CorrectDistortion,1);
    end
    
    if MV.P.Mirror1; cData = flipud(cData); end % Corresponds to X if camera on its base
    if MV.P.Mirror2; cData = fliplr(cData); end % Corresponds to Y if camera on its base
    if MV.P.Transpose; cData = cData'; end
    
    % OUTPUT FRAME TO FRAME DIFFERENCE IN BRIGHTNESS
    Res = sort(cData(:));
    CorrFact = mean(Res(end-round(0.01*length(Res)):end));
    CorrFact = CorrFact/220;
    cData = double(cData)/CorrFact;
    
    %% SHOW OTHER INFORMATION (WHISKER FITS / CONTACTS)
    try delete(MV.LastPH); end
    MV.LastPH = [];
    if ~isnan(cFrame)
      if ~isempty(MV.CurrentFrame) & ~isempty(MV.Annotations{MV.CurrentFrame}) ...
          & (isstruct(MV.Annotations{MV.CurrentFrame}) | ~isnan(MV.Annotations{MV.CurrentFrame}))
        cCurves = MV.Annotations{MV.CurrentFrame}.Clicks;
        if ~isempty(cCurves)
          Counts = [0,0,0];
          for i=1:length(cCurves)
            switch cCurves(i).ClickType
              case 'left'; cColor = [1,0.5,0.5]; cInd = 1; % Red chosen to mark female
              case 'right'; cColor = [0.5,0.5,1];  cInd = 2; % Blue chosen to mark male
              case 'center'; cColor = [0.5,1,0.5];  cInd = 3; % Third animal in triad, can either be male or female
              otherwise cColor = [1,1,1];
            end
            Counts(cInd) = Counts(cInd)+1;
            if Counts(cInd)>length(MV.P.Markers) Counts(cInd) = length(MV.P.Markers); end
            MV.LastPH(end+1) = plot(MV.AH.Video,...
              cCurves(i).X,cCurves(i).Y,MV.P.Markers{Counts(cInd)},...
              'Color',cColor,'MarkerSize',6,'LineWidth',2,'Hittest','off');
          end
        end
      end
    end
  else
    try delete(MV.LastPH); end
    
    cData = zeros(R.(MV.VideoName).Dims);
    cFrame = NaN;
  end

  set(MV.GUI.Video.FrameH,'CData',cData'); % NOTE: Transpose, due to keeping image in XY before plotting
  set(MV.GUI.Video.TitleH,'String',['Video : Trial ',num2str(Trial),' | Frame : ' num2str(cFrame),' '])
 
function MV_showLocalization % Show Localization Results from Cam64
  global MV R;
  
  Trial = find(MV.CurrentTime>MV.LocalizationStarts,1,'last');
  if ~isempty(Trial) && Trial > 0
    [cDiff,cFrame] = min(abs(R.Localization.Data(Trial).Data.Time(:,1) - MV.CurrentTime));
    
    if ~isempty(cFrame) && cDiff < 0.02 % if within 50ms of an existing frame
      cData = squeeze(R.Localization.Data(Trial).Data.Frames(:,:,cFrame)); 
      cMedian = median(cData(:)); cMax = max(cData(:));
      cData = cData-cMedian; cData = cData/(0.3*(cMax-cMedian));
      if max(cData(:))>1 cData = 1.2*cData/(max(cData(:))); end

      cX = R.Localization.Data.Data.Location(cFrame,1);
      cY = R.Localization.Data.Data.Location(cFrame,2);
    else 
      cData = zeros(R.Localization.Dims);  cX = NaN; cY = NaN;
    end
  else
    cData = zeros(R.Localization.Dims);  cX = NaN; cY = NaN;    
  end
  set(MV.GUI.Localization.MaxH,'XData',cX,'YData',cY);
  
  %   if ~isempty(Trial) & isfield(R.Localization.Data(Trial).Data,'Vec')
  %     cVec = R.Localization.Data(Trial).Data.Vec(cFrame,:);
  %     cVec = 0.01*cVec;
  %     set(MV.GUI.Localization.DirH,'XData',cX + [-cVec(1),cVec(1)] ,'YData',cY+ [-cVec(2),cVec(2)]);
  %   end
  cImage = cat(3,cData,zeros([size(cData),2])); %
  set(MV.GUI.Localization.FrameH,'CData',cImage,'AlphaData',(cImage(:,:,1).^2)/2);
  
function MV_showAVLocalization
global MV R

for i=1:length(MV.P.Vocs) 
  Locs(i,:) = MV.P.Vocs(i).Location;
  TimesVoc(i) = (MV.P.Vocs(i).Start + MV.P.Vocs(i).Stop)/2;
  [~,Frames(i)] = min(abs(R.PointGrey.Data.Data.Time(:,1)-TimesVoc(i)));
  C{i}=MV.Annotations{Frames(i)}.Clicks;
  cC=C{i}; 
  CTs = {cC.ClickType}; 
  Ts = {cC.Type}; 
  cInd = strcmp(CTs,'right') & strcmp(Ts,'Snout'); 
  MalePos(i,:) = [cC(cInd).X,cC(cInd).Y]; 
end

figure(1000); clf
[~,AH]=axesDivide(2,1,'c');
for iA =1:length(AH)
  axes(AH(iA)); hold on;
  plot(Locs(:,iA),MalePos(:,iA),'.b');
  plot([-0.2,0.2],[-0.2,0.2],'Color',[0.5,0.5,0.5])
  axis([-0.2,0.2,-0.2,0.2]);
  axis square; grid on; 
end
MEA = median(abs(sqrt(sum((Locs - MalePos).^2,2))));
fprintf(['MEA = ',num2str(MEA,4),'m\n']);
keyboard

function Frame = MV_convertTime2Frame(Time)
  global MV R
  Trial = find(Time>MV.VideoStarts,1,'last');
  Frame = find(R.(MV.VideoName).Data(Trial).Data.Time(:,1) >= Time,1,'first');
  
function Time = MV_convertFrame2Time(Frame)
  global MV R
  LastStart = 0; 
  Time = NaN;
  for iT=1:length(R.(MV.VideoName).Data)
    cFrames = LastStart + size(R.(MV.VideoName).Data(iT).Data.Time,1);
    if Frame<=cFrames; 
      Time  = R.(MV.VideoName).Data(iT).Data.Time(Frame+LastStart,1);
      break;
    else 
      LastStart = LastStart + cFrames;
    end
  end
  
function MV_gotoFrame
    global MV R
    cFrame = input('Enter frame: ');
    cTime = MV_convertFrame2Time(cFrame);
    MV_showData([],[],'set',cTime)
  
function MV_loadAnnotations(O,E,Option)
  global MV R
  if nargin<3 Option = 'Ask'; end
  cAnimal = R.General.Parameters.General.Animal;
  cRecording = R.General.Parameters.General.Recording;
  [~,DD] = C_getDir('Animal',cAnimal,'Recording',cRecording);
  DefaultPath = [fullfile(DD.DB,'Results','Video'),filesep];
  FallbackPath = ['C:\Users\Lab\Desktop\'];
  switch Option
    case 'Manual'
      File = [cAnimal,'_R',num2str(cRecording),'_vocvideotracking.mat'];
      if exist([DefaultPath,File],'file'); BasePath = DefaultPath; else clear File; BasePath = FallbackPath; end
    case 'Max'
      File = [cAnimal,'_R',num2str(cRecording),'_videotracking_max_corrected.h5.mat'];
      if exist([DefaultPath,File],'file'); BasePath = DefaultPath; else clear File; BasePath = FallbackPath; end
    case 'Automatic'
      File = [cAnimal,'_R',num2str(cRecording),'_videotracking.h5*'];
      Files = dir([DefaultPath,File]);
      if ~isempty(Files); BasePath = DefaultPath; else clear File; BasePath = FallbackPath; end
    otherwise
      if exist(DefaultPath,'Dir') BasePath = DefaultPath; else BasePath = FallbackPath; end
  end
  if ~exist('File','var');  [File,BasePath] = uigetfile([BasePath,'*']); end
  if ~File return; end
  Pos = find(File=='.');
  Ext = File(Pos(1)+1:end);
  switch Ext;
    case 'mat'; % Load in native format
      tmp = load([BasePath,File]); FN = fieldnames(tmp);
      if sum(strcmp(FN,'Annotations')) % New Format
        MV.Annotations = tmp.Annotations;
      else % Old Format
        CBF = tmp.CurvesByFrame;
        Ind = find(~cellfun(@isempty,CBF));
        A = cell(size(CBF));
        for iF=1:length(Ind)
          cFrame = Ind(iF);
          cCBF = CBF{cFrame};
          for iC=1:length(cCBF)
            cXPoint = cCBF{iC}{1};
            cYPoint = cCBF{iC}{2};
            cClickType = cCBF{iC}{4};
            cTime = MV_convertFrame2Time(cFrame);
            switch iC
              case {1,3}; cType = 'Snout';
              case {2,4}; cType = 'HeadCenter';
            end
            [~,cXPixel] = min(abs(MV.SetupInfo.Video.XPos - cXPoint));
            [~,cYPixel] = min(abs(MV.SetupInfo.Video.YPos - cYPoint));
            A{cFrame}.Clicks(iC) = ...
              struct('X',cXPoint,'Y',cYPoint,'XPixel',cXPixel,'YPixel',cYPixel,...
              'Time',cTime,'Frame',cFrame,...
              'Type',cType,'ClickType',cClickType);
          end
        end
        %warning('Shifting Annotations by hand to correct for earlier alignment error.');
        MV.Annotations = A;
      end
      
    case 'csv';
      [Items,Data] = loadCSV([BasePath,File],'Separator',',');
      for iA=1:size(Items,1)
        cFrame = str2num(Items{iA,1}(4:end-4));
        cAnno = MV.Annotations{cFrame};
        if isstruct(cAnno) cClickInd = length(cAnno.Clicks) + 1; else cClickInd = 1; end
        cTime = MV_convertFrame2Time(cFrame);
        cTypeRaw = Items{iA,3};
        Pos = find(cTypeRaw=='_');
        cType = cTypeRaw(1:Pos-1);
        cClickType = cTypeRaw(Pos+1:end);
        cXPixel = Items{iA,4};
        cYPixel = Items{iA,5};
        cXPoint = MV.SetupInfo.Video.XPos(cXPixel);
        cYPoint = MV.SetupInfo.Video.YPos(cYPixel+1);
        MV.Annotations{cFrame}.Clicks(cClickInd) = ...
          struct('X',cXPoint,'Y',cYPoint,'XPixel',cXPixel,'YPixel',cYPixel,...
          'Time',cTime,'Frame',cFrame,...
          'Type',cType,'ClickType',cClickType);
      end
      
    case 'h5'; % DeepLabCut Format
      [Frames,Nums,Tags,Types] = C_loadAnnotationsDLC([BasePath,File]);
      NFields = length(Types);
      TypesInd = 1:NFields;
      
      % ASSIGN DATA TO INTERNAL FORMAT
      for iF=1:length(Frames)
        cFrame = Frames(iF);
        cFrames = mat2cell(repmat(cFrame,1,length(Types)),1,ones(1,length(Types)));
        if ~mod(cFrame,2000) fprintf([num2str(cFrame),' ']); end
        cTime = MV_convertFrame2Time(cFrame);
        if ~isnan(cTime)
          cTimes = mat2cell(repmat(cTime,1,length(Types)),1,ones(1,length(Types)));
          cXPixelsTmp = Nums(cFrame,(TypesInd-1)*3+2);
          cYPixelsTmp = 512-Nums(cFrame,(TypesInd-1)*3+1);
          cXPixels = mat2cell(cXPixelsTmp,1,ones(1,length(Types)));
          cYPixels = mat2cell(cYPixelsTmp,1,ones(1,length(Types)));
          cXPos = mat2cell(MV.SetupInfo.Video.Pix2Pos(cXPixelsTmp,'X'),1,ones(1,length(Types)));
          cYPos = mat2cell(MV.SetupInfo.Video.Pix2Pos(cYPixelsTmp,'Y'),1,ones(1,length(Types)));
          %cTypeShifts = mat2cell(,1,ones(1,length(Types)));
          MV.Annotations{cFrame}.Clicks = ...
            struct('X',cXPos,'Y',cYPos,'XPixel',cXPixels,'YPixel',cYPixels,...
            'Time',cTimes,'Frame',cFrames,...
            'Type',Types,'TypeShift',[],'VocLocation',NaN,'ClickType',Tags);
        end
      end
      fprintf('... Done\n');
      
    case 'h5.mat' % Converted DLC Format
      tmp = load([BasePath,File]); D = tmp.cData;
      % ASSIGN DATA TO INTERNAL FORMAT
      FN = setdiff(fieldnames(D),{'Times','Frames'});
      for iF = 1:length(FN)
        Pos = find(FN{iF}=='_');
        Types{iF} = FN{iF}(1:Pos-1); Tags{iF} = FN{iF}(Pos+1:end);
      end
      NFields = length(FN);
      for iF=1:length(D.Frames)
        cFrame = D.Frames(iF);
        cFrames = mat2cell(repmat(cFrame,1,NFields),1,ones(1,NFields));
        if ~mod(cFrame,2000) fprintf([num2str(cFrame),' ']); end
        cTime = MV_convertFrame2Time(cFrame);
        if ~isnan(cTime)
          cTimes = mat2cell(repmat(cTime,1,NFields),1,ones(1,NFields));
          DataCenter = MV.SetupInfo.Video.Offset{1}(1:2); Zoom = 1; Strength =0.3;
          for iF = 1:length(FN)
            %             Data = D.(FN{iF})(cFrame,1:2)';
            %             DataCorrected = C_correctImageDistortion(Data,Strength,Zoom,'Markers',DataCenter);
            %             cXPos{iF} = DataCorrected(1);
            %             cYPos{iF} = DataCorrected(2);
            
            cXPos{iF} = D.(FN{iF})(cFrame,1);
            cYPos{iF} = D.(FN{iF})(cFrame,2);
          end
          cXPixels = mat2cell(repmat(NaN,1,NFields),1,ones(1,NFields));
          cYPixels = mat2cell(repmat(NaN,1,NFields),1,ones(1,NFields));
          MV.Annotations{cFrame}.Clicks = ...
            struct('X',cXPos,'Y',cYPos,'XPixel',cXPixels,'YPixel',cYPixels,...
            'Time',cTimes,'Frame',cFrames,...
            'Type',Types,'TypeShift',[],'VocLocation',NaN,'ClickType',Tags);
        end
      end
      fprintf('... Done\n');      
      
    otherwise error('File format not known');
  end
  disp(['Loading previously collected Annotations from [ ',BasePath,File,' ]']);
  MV_showData([],[],'redraw');

function MV_saveAnnotations(O,E)
  global MV R
  if ~isfield(MV,'Annotations') return; end
  Res.Annotations = MV.Annotations;
  AnnotationsExist = sum(~cellfun(@isempty,Res.Annotations));
  
  if AnnotationsExist
    Format = input('Choose Format: [M]atlab or [C]SV: ','s');
    Format = lower(Format);
    switch Format; case 'c'; Extension = 'csv'; case 'm'; Extension = 'mat'; otherwise return; end
    MV.FileNameAnnotations = input('Enter Filename for Annotations [relative to current Path, no extension] : ','s');
    if isempty(MV.FileNameAnnotations) return; end
    MV.FormatAnnotations = Extension;
    MV.FileNameAnnotations = [pwd,filesep,MV.FileNameAnnotations,'.',Extension];
    disp(['Saving to ',MV.FileNameAnnotations,'.']);
    
    switch MV.FormatAnnotations
      case 'mat'; % MATLAB FILE
        save(MV.FileNameAnnotations,'-struct','Res');
      case 'csv'; % CSV FILE
        Frames = MV_getTrackedFrames; 
        NTrackPoints = length(MV.Annotations{Frames(1)}.Clicks);
        XY = zeros(length(Frames),NTrackPoints*2);
        FID = fopen(MV.FileNameAnnotations,'w');
        for iF=1:length(Frames)
          cNTrackPoints = length(MV.Annotations{Frames(iF)}.Clicks);
          if cNTrackPoints ~=NTrackPoints; 
            fprintf(['In Frame ',num2str(Frames(iF)),' : Actual # of Trackpoints: ',num2str(cNTrackPoints),', expected ',num2str(NTrackPoints),'\n']);
          end
          for iT=1:NTrackPoints
            cClick = MV.Annotations{Frames(iF)}.Clicks(iT);
            cType = [cClick.Type,'_',cClick.ClickType];
            if isfield(cClick,'TypeShift'); cType = [cType,'_',cClick.TypeShift]; end
            fprintf(FID,'%s,%s,%s,%d,%d\n',['img',sprintf('%06d',Frames(iF)),'.jpg'],'default',cType,cClick.XPixel,513-cClick.YPixel);
          end
        end
        fclose(FID);
    end
  end
  
function Frames = MV_getTrackedFrames
  global MV;
  Frames = find(~cellfun(@isempty,MV.Annotations));
  RemFrames = [];
  for iF=1:length(Frames)
    cFrame = Frames(iF); cAnnotation = MV.Annotations{cFrame};
    if ~isstruct(cAnnotation) & isnan(cAnnotation); RemFrames(end+1) = cFrame; end
  end
  Frames = setdiff(Frames,RemFrames);
  
function MV_scrobblor(O,E)

  global MV R;
  
  AxPos = get(MV.AH.Scrobble,'Position');
  FigPos = get(MV.FIG,'Position');
  Pixels = AxPos(3)*FigPos(3);
  TimeTotal = MV.StopTime-MV.StartTime;
  TimePerPixel = TimeTotal/Pixels;
  AxisStartPixels =  FigPos(1) + AxPos(1)*FigPos(3);
  % LastPoint = [-inf,-inf];
  
  SelType = get(MV.FIG,'SelectionType');
  switch SelType
    case 'normal';
      global Scrobbling_ ; Scrobbling_ =1;
      MV.P.ComputeLocation = 0;
      while Scrobbling_
        CurrentPoint = get(0,'PointerLocation');
        CurrentTime  = TimePerPixel*(CurrentPoint(1)-AxisStartPixels) + MV.StartTime;
        MV_showData([],[],'set',CurrentTime);
        global Scrobbling_ ;
        pause(0.2);
      end
      
    case {'alt','extend'};
      MV.P.ComputeLocation = 1;
      global Scrobbling_ ; Scrobbling_ =0;
  end
  
function MV_setTime(O,E,Time)
  global MV R;
  
  if nargin==3 
    NewTime = Time;
  else
    Position = get(gca,'CurrentPoint');
    NewTime = MV.CurrentTime + Position(1);
    
    SelType = get(MV.FIG,'SelectionType');
    switch SelType
      case {'alt','extend'};
        MV.P.ComputeLocation = 1;
      otherwise
        MV.P.ComputeLocation = 0;
    end;
  end
  
  MV_showData([],[],'set',NewTime);
 
function [AnimalPositionsVsTime,AnimalTime,AnimalPositions] = MV_computeAnimalPosition
    global MV R;
    Pos = C_estimateAnimalPosition('Data',R,'BetweenSensors',0,'SR',10);
    AnimalPositionsVsTime = Pos.AnimalPositions;
    AnimalTime = Pos.Time;
    AnimalPositions = Pos.Sensors;
    NPos =size(AnimalPositionsVsTime,1);
    AnimalPositions = linspace(1/(2*NPos), 1-1/(2*NPos),NPos);    
    
function MV_setCLimit(O,E,Opt,CLimits)
global MV R;
if ~exist('Opt','var') Opt = ''; end

switch Opt
  case 'set'; MV.CLimits = CLimits;
  case 'setedit'; 
    if O==MV.GUI.CLimMin LimInd = 1; else LimInd = 2; end 
    MV.CLimits(LimInd) = str2num(get(O,'String'));
  case 'setvisual';
    NewLimit = get(MV.AH.Colorbar,'CurrentPoint');
    NewLimit = round(255*NewLimit(1,2));
    SelType = get(MV.FIG,'SelectionType');
    switch SelType
      case 'normal'; MV.CLimits(1) = NewLimit;
      case 'alt'; MV.CLimits(2) = NewLimit;
    end
  otherwise % for empty case
end

MV.CLimits = sort(MV.CLimits);
CY = linspace(MV.P.CMin,MV.P.CMax,MV.NColSteps);
NBelow =  round(MV.NColSteps*(MV.CLimits(1)-MV.P.CMin)/(MV.P.CMax-MV.P.CMin));
NAbove = round(MV.NColSteps*(MV.P.CMax-MV.CLimits(2))/(MV.P.CMax-MV.P.CMin));
CY = [zeros(1,NBelow) , linspace(0,1,MV.NColSteps-NBelow-NAbove) , ones(1,NAbove)]' ;
set(MV.GUI.CIH,'CData',CY,'YData',CY,'XData',1); caxis(MV.AH.Colorbar,[0,1]);
set(MV.AH.Colorbar,'YLim',CY([1,end]),'YTick',[]);
if MV.CLimits(1) == MV.CLimits(2) MV.CLimits = MV.CLimits+[0,1]; end
caxis(MV.AH.Video,MV.CLimits);

function MV_DetectVocs
    global MV R;
    if ~isfield(R.AnalogIn,'FromWAV')
      cAnimal = R.General.Parameters.General.Animal;
      cRecording = R.General.Parameters.General.Recording;
      MV.P.Vocs = VocCollector('Animals',{cAnimal},'Recording',cRecording);
    else
      MV.P.Vocs = VocCollector('DataSource','Direct','Data',R.AnalogIn.Data.Data.Analog,'SRAI',R.AnalogIn.SRAI);
    end
