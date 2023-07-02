%% Sample script for demonstrating the use of HyVL on a sample data set that is part of the repository

% ADD THE REQUIRED PATH, INCLUDING THE HELPER FILES
addpath(genpath('.');

% LOAD DATA & SETUP INFORMATION
load('SampleData.mat');

% DETECT VOCALIZATIONS
D = R.AnalogIn.Data.Data.Analog; 
SRUSM4 = 1/diff(R.AnalogIn.Data.Data.Time(1:2));
Vocs = VocCollector('Data',D,'SRAI',SRUSM4);

% LOCALIZE VOCALIZATIONS USING USM4 DATAS
MP = Setup.Audio.MicrophonePositions;
SH = Setup.SourceHeight;
Vocs = VocAnalyzer(Vocs,'MicrophonePositions',MP,'SourceHeight',SH);

% LOCALIZE VOCALIZATIONS USING CAM64 DATA
D = R.Cam64.Data.Data.Analog;
SRCam64 = 1/diff(R.Cam64.Data.Data.Time(1:2));
Localization = VocLocalizerCam64('Data',D,'SR',SRCam64,'Vocs',Vocs,'Setup',Setup);
R.Localization =  C_convertCam64ResultsToController('Data',Localization.Vocs);

% SHOW DATA, TRACKING AND LOCALIZATIONS
global R;
MultiViewer('PhysicalSetup',Setup,'Vocs',Vocs,'LoadAnnotations',0,'CorrectDistortion',0.34,'ShowSensors',0);

%% %%%%%%%%%%%%%% WAS ONLY USED FOR EXPORTING DATA FOR ARCHIVE %%%%%%%%%%%%%
function copyToRepo

SourcePath = '~/GoogleDrive/Code/Controller/';
TargetPath = '~/GoogleDrive/Projects/ART_USVLocCam64/github/HyVL/';
Files = {...
  'Access/Vocalizations/Cam64/BasicHyVLPipeline.m','';...
  'Access/Vocalizations/VocLocalizer.m','Main';...
  'Access/Vocalizations/Cam64/VocLocalizerCam64.m','Main';...
  'Access/Vocalizations/VocAnalyzer.m','Main';...
  'Access/Vocalizations/Cam64/C_convertCam64ResultsToController.m','Main';...
  'Plot/Helpers/axesDivide.m','Helpers';...
  'Helpers/parsePairs.m','Helpers';...
  'Helpers/checkField.m','Helpers';...
  'Access/MultiViewer/MultiViewer.m','Main',...
};

for iF=1:length(Files)
  copyfile([SourcePath,Files{iF,1}],[TargetPath,Files{iF,2}]);
end



function reduceTimeRange

RR = C_loadRecording('Animal','mouse82','Recording',45,'VideoLoadMode','Full','Cam64LoadMode','Full');

TimeRange = [38.6,48.6];

Modules = {'NIDAQ','AnalogIn','Cam64','PointGrey'};
R.General = RR.General;
% REDUCE RANGE IN MODULES
for iM = 1:length(Modules)
  cM = Modules{iM};
  R.(cM) = RR.(cM);
  Inds = find(RR.(cM).Data.Data.Time>TimeRange(1) & RR.(cM).Data.Data.Time<TimeRange(2));
  switch cM
    case 'PointGrey';  R.(cM).Data.Data.Frames = RR.(cM).Data.Data.Frames(:,:,:,Inds);
    case 'Cam64'; R.(cM).Data.Data.Analog = RR.(cM).Data.Data.Analog{1}(Inds,:);
    otherwise R.(cM).Data.Data.Analog = RR.(cM).Data.Data.Analog(Inds,:);
  end
  R.(cM).Data.Data.Time = RR.(cM).Data.Data.Time(Inds,:) - RR.(cM).Data.Data.Time(Inds(1),:);
end
R.PointGrey.Data.Data.Time = R.PointGrey.Data.Data.Time+ 0.2; % Compensate dropped frames for this particular video
R.NIDAQ.Data.CameraStartTime(1) = 0;  
R.NIDAQ.Data.LocalizationStartTimes = 0;
R.PointGrey.Data.Sync.SyncTimeRel = 0;

Setup = C_getSetupInfo('PhysicalSetup','Platform2DCam64Flea3');

% SAVE FILE
save('SampleData.mat','R','Setup');