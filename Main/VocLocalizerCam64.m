function RR = VocLocalizerCam64(varargin)
% Sample for Methods figure: M84_R64_F12694. 

P = parsePairs(varargin);
checkField(P,'Data')
checkField(P,'SR');
checkField(P,'Setup')
checkField(P,'Vocs')
checkField(P,'FRangeMax',50000)
checkField(P,'FMax',90000) % Limit the range of frequencies include, since noise level grows with Frequency
checkField(P,'PhysicalSetup','Platform2DCam64Flea3')
checkField(P,'PeakMethod','Average')
checkField(P,'Verbose',0);
checkField(P,'FIG',1);
checkField(P);

XStep.Coarse = 0.01; YStep = XStep; P.Simulation = 0;
switch gpuDeviceCount
  case 0; Method = 'DelaySum';
  otherwise Method = 'DelaySumGPU';
end
LocalOpts = {'Method',Method,'TRange',[0],'Verbose',P.Verbose};

% SETUP GENERAL PARAMETERS
SR.Cam64 = P.SR;
Time.Cam64 = [0:size(P.Data,1)-1]/SR.Cam64; 
SI = P.Setup;
XRange.Coarse = [-0.25,0.25]; YRange.Coarse = [-0.2,0.2];
YStep.Coarse = XStep.Coarse; XStep.Fine = 0.001; YStep.Fine = 0.001;
Resolutions = {'Coarse','Fine'};
RelRangeFine = [-0.015,0.015];
NFFT = 256;

Mics = C_mapCam64;  NMics = size(Mics.LocByMic,1);
ZDist.Cam64 = SI.Cam64.MicrophonePosition{1}(3) - SI.SourceHeight;
MicPos.Cam64 = [... % Note that X and Y are switched here, since Cam64 rotated in the setup
  Mics.LocByMic(:,2) - Mics.Center(2) + SI.Cam64.MicrophonePosition{1}(1), ... % X
  Mics.LocByMic(:,1) - Mics.Center(1) + SI.Cam64.MicrophonePosition{1}(2), ... % Y
  zeros(NMics,1)];

% PREPARE FOR SAVING
RR = [];

% PREPARE FIGURE FOR PLOTTING
if P.FIG;  figure(P.FIG); clf; [~,AH] = axesDivide(3,1,'c'); end

%% ==============================================================
for iV = 1 : length(P.Vocs)
  printupdate([num2str(iV),'/',num2str(length(P.Vocs))],iV==1);
  cVoc = P.Vocs(iV);

  % SET MATCHED TIME AND FREQUENCY RANGES
  cTRange = [cVoc.Start,cVoc.Stop];
  cFRange.Voc = [cVoc.FMin-1000,cVoc.FMax+1000];
  cFRange.Voc = [40000,90000];
  cFRange.Cam64 = cFRange.Voc;
  if cFRange.Cam64(1) < P.FMax
    if cFRange.Cam64(2) > P.FMax;  cFRange.Cam64(2) = P.FMax + 5000;  end
    if cFRange.Cam64(1) > cFRange.Cam64(2); cFRange.Cam64(1) = cFRange.Cam64(2) - 5000; end
  end
  if diff(cFRange.Cam64)>P.FRangeMax; cFRange.Cam64 = cVoc.FMean + P.FRangeMax*[-0.5,0.5]; end
  cFRange.Cam64(2) = min(floor(SR.Cam64/2)-1,cFRange.Cam64(2));

  cIndAnalog = find(Time.Cam64>=cVoc.Start & Time.Cam64<cVoc.Stop) ;
  cData.Cam64 = double(P.Data(cIndAnalog(1):cIndAnalog(end),:));
  cData.Cam64 = cData.Cam64./std(cData.Cam64);
  
  cDevice = 'Cam64'; clear S; cAH = AH(1,1); cla(cAH); hold(cAH, 'on');
  if ~isempty(cData.(cDevice))
    cTime = (cVoc.Time-cVoc.Time(1))*1000;
    for iC=1:size(cData.(cDevice),2)
      [S(:,:,iC),F,T] = HF_specgram(cData.(cDevice)(:,iC),NFFT,SR.(cDevice),[10000,125000],NFFT/2,0,1);
    end
    S = mean(S,3); T = T*1000; F = F/1000;
    h=pcolor(cAH,T,F,S); set(h,'EdgeC','None');
    plot(cAH,cTime,repmat(cVoc.FMin/1000,size(cTime)),'g--');
    plot(cAH,cTime,repmat(cFRange.(cDevice)(1)/1000,size(cTime)),'g-');
    plot(cAH,cTime,repmat(cVoc.FMax/1000,size(cTime)),'r--');
    plot(cAH,cTime,repmat(cFRange.(cDevice)(2)/1000,size(cTime)),'r-');
    xlabel(cAH,'Time [ms]'); title(cAH,cDevice);
    ylim(cAH,[0,125]);% xlim(cAH, cTime([1,end]))
    ylabel(cAH,'Freq. [kHz]');
    text(-1.15,1,['TStart=',num2str(cVoc.Start)],'Units','n','FontSize',7);
    text(-1.15,0.9,['TStop=',num2str(cVoc.Stop)],'Units','n','FontSize',7);
  end
  drawnow;


  % LOOP OVER ESTIMATION LEVELS
  clear cRes
  for iR = 1:length(Resolutions)
    cResolution = Resolutions{iR};
    cDevice = 'Cam64';
    if ~isempty(cData.(cDevice))
      B = C_BeamformingTime('Data',cData.(cDevice),'SR',SR.(cDevice),...
        'MicrophonePositions',MicPos.(cDevice),...
        'XRange',XRange.(cResolution),'YRange',YRange.(cResolution),...
        'XStep',XStep.(cResolution),'YStep',YStep.(cResolution),'ZDist',ZDist.(cDevice),...
        'TBin',size(cData.(cDevice),1)/SR.(cDevice),'FRange',cFRange.(cDevice),LocalOpts{:});
      cRes.(cResolution).Time = cVoc.Start;
    end

    % Normalize Beamforming Results
    switch cResolution
      case 'Coarse'
        Median.(cDevice) =  median(B.ProjectionXY(:));
        Std.(cDevice) = std(B.ProjectionXY(:) - Median.(cDevice));
    end
    B.ProjectionXY = B.ProjectionXY - Median.(cDevice);
    B.ProjectionXY = B.ProjectionXY/Std.(cDevice);
    cRes.(cResolution).(cDevice) = B;

    cRes.(cResolution).(cDevice).Location = LF_getPeak(cRes.(cResolution).(cDevice),cResolution,P.PeakMethod);
    cRes.(cResolution).(cDevice).SNR = max(B.ProjectionXY(:))/std(B.ProjectionXY(:));

    % FUSE ESTIMATES ACROSS DEVICES
    cRes.(cResolution).Combined = cRes.(cResolution).Cam64;
    
    % DETERMINE LOCATION
    cRes.(cResolution).Combined.Location = LF_getPeak(cRes.(cResolution).Combined,cResolution,P.PeakMethod);

    if P.FIG % SHOW BEAMFORMING RESULTS
      CM = HF_colormap({[0,0,1],[1,1,1],[1,0,0]},[-1,0,1]);
      if ~isempty(cData.(cDevice))
        iA = iR + 1; cAH = AH(iA); cla(cAH); %axes(cAH);
        imagesc(cAH,cRes.(cResolution).(cDevice).x,cRes.(cResolution).(cDevice).y,cRes.(cResolution).(cDevice).ProjectionXY);
        set(cAH,'YD','n'); title(AH(iA),[cDevice,' ',cResolution]);  hold(cAH,'on');
        plot(cAH,cRes.(cResolution).(cDevice).Location(1),cRes.(cResolution).(cDevice).Location(2),'go','Markersize',10);
        axis(cAH,'tight'); colorbar(cAH);  colormap(cAH,CM);
      end
      drawnow;
    end

    % DETERMINE COARSE MAXIMUM AS STARTING POINT FOR FINE ESTIMATE
    if strcmp(cResolution,'Coarse')
      XRange.Fine = cRes.Coarse.Combined.Location(1) + RelRangeFine ;
      YRange.Fine = cRes.Coarse.Combined.Location(2) + RelRangeFine;
    end

  end % RESOLUTION
  RR.Vocs(iV).Start = cVoc.Start;
  RR.Vocs(iV).Stop = cVoc.Stop;
  RR.Vocs(iV).ResCoarse = cRes.Coarse.Combined;
  RR.Vocs(iV).ResFine = cRes.Fine.Combined;

end % VOCALIZATIONS


function PosXY = LF_getPeak(cRes,cResolution,PeakMethod);

switch cResolution
  case 'Coarse'; cPeakMethod = 'Max';
  case 'Fine'; cPeakMethod = PeakMethod;
end

switch cPeakMethod
  case 'Max'
    [i1,i2] = find(cRes.ProjectionXY == max(cRes.ProjectionXY(:)));
    PosXY = [cRes.x(i2),cRes.y(i1)];
  case 'Average'
    cP = cRes.ProjectionXY; cP(cP<0) = 0; cP = cP.^2; cP = cP./sum(cP,'all');
    [cX,cY] = meshgrid(cRes.x,cRes.y);
    PosXY = [sum(cP.*cX,'all'),sum(cP.*cY,'all')];
end