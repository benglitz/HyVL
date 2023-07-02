function [Localization,R] = C_convertCam64ResultsToController(varargin)
% Strategy:
% - Output: Images compatible with MultiViewer and Absolute locations
% TODO:
% - Integrate scale/shift if necessary
% - Adjust to minimize error on recordings with single mouse
% - Saving the results to the Vocalization file, so that SC_VocAnalyzer and C_analyzeSnout.. can pick them up.

P = parsePairs(varargin);
checkField(P,'Animal',[])
checkField(P,'Recording',[])
checkField(P,'Backend','Local')
checkField(P,'Data',[])
checkField(P,'Devices',{'Cam64'})
checkField(P,'PeakMethod','Average')
checkField(P,'ComputeDirection',0)
checkField(P,'Setup','Platform2DCam64Flea3')
checkField(P,'SetupInfo',[])
checkField(P,'FIG',0)
checkField(P)

% LOAD DATA IF NOT PASSED FROM OUTSIDE
if isempty(P.Data)  I = C_getRecInfo('Animal',P.Animal,'Recording',P.Recording);
  cResultsFile = ['Localization_',P.Devices{:},'_',P.Backend,'.mat'];
  cResultsFileLocal = [I.Local.Results.Path,'Vocalizations',filesep,cResultsFile];
  R = load(cResultsFileLocal); Vocs = R.RR.Vocs;
else 
  Vocs = P.Data;
end

% COLLECT SETUP INFORMATION IF NOT PASSED FROM OUTSIDE
if isempty(P.SetupInfo)
  SI = C_getSetupInfo('PhysicalSetup',P.Setup);
else
  SI = P.SetupInfo;
end

% SETUP FOR BEAMFORMING
[XCoarse,YCoarse] = meshgrid(Vocs(1).ResCoarse.x,Vocs(1).ResCoarse.y);
XResCoarse = diff(Vocs(1).ResCoarse.x(1:2));
YResCoarse = diff(Vocs(1).ResCoarse.y(1:2));
XRefinedLimits = [Vocs(1).ResCoarse.x(1)-XResCoarse/2 , Vocs(1).ResCoarse.x(end)+XResCoarse/2]; 
YRefinedLimits = [Vocs(1).ResCoarse.y(1)-YResCoarse/2 , Vocs(1).ResCoarse.y(end)+YResCoarse/2]; 
XResFine = diff(Vocs(1).ResFine.x(1:2));
YResFine =  diff(Vocs(1).ResFine.y(1:2));
[XRefined,YRefined] = meshgrid(...
  [XRefinedLimits(1):XResFine:XRefinedLimits(2)],...
  [YRefinedLimits(1):YResFine:YRefinedLimits(2)]);
Localization.X = XRefined(1,:);
Localization.Y = YRefined(:,1);
Localization.FromCam64 = 1;

% LOOP OVER VOCALIZATIONS
for iV=1:length(Vocs)
  printupdate([num2str(iV),'/',num2str(length(Vocs))],iV==1);
  % FUSE COARSE AND FINE (Oversample and combine)
  if isempty(Vocs(iV).ResCoarse) continue; end
  cProjectionRefined = interp2(XCoarse,YCoarse,Vocs(iV).ResCoarse.ProjectionXY,XRefined,YRefined,'makima');
  cProjectionFine = Vocs(iV).ResFine.ProjectionXY;
  XFine = Vocs(iV).ResFine.x;
  YFine = Vocs(iV).ResFine.y;
  [~,XIndStart] = min(abs(XRefined(1,:)-XFine(1)));
  [~,XIndEnd] = min(abs(XRefined(1,:)-XFine(end)));
  XInd = [XIndStart : XIndEnd];
  [~,YIndStart] = min(abs(YRefined(:,1)-YFine(1)));
  [~,YIndEnd] = min(abs(YRefined(:,1)-YFine(end)));
  YInd = [YIndStart : YIndEnd];
  if length(YInd)~=size(cProjectionFine,1)
    if YIndStart == 1
      YIndFine = [size(cProjectionFine,1)-length(YInd)+1:size(cProjectionFine,1)];
    else
      YIndFine = [1:length(YInd)];
    end
  else
    YIndFine = [1:size(cProjectionFine,1)];
  end
  if length(XInd)~=size(cProjectionFine,2)
    if XIndStart == 1
      XIndFine = [size(cProjectionFine,2)-length(XInd)+1:size(cProjectionFine,2)];
    else
      XIndFine = [1:length(XInd)];
    end
  else
    XIndFine = [1:size(cProjectionFine,2)];
  end
  cProjectionRefined(YInd,XInd) = cProjectionFine(YIndFine,XIndFine);
  
  cR.Frames(:,:,iV) = cProjectionRefined;
  cR.Time(iV,1) = (Vocs(iV).Start + Vocs(iV).Stop)/2 ;
  
  
  % 
  MAX = max(cProjectionFine(:));
  MIN = min(cProjectionFine(:));
  if MAX > MIN
    [i1,i2] = find(cProjectionFine ==MAX);
    Center = round([mean(i2);mean(i1)]); % [X;Y]
    
    switch P.PeakMethod
      case 'Max' % Just take the point of maximal value
        cR.Location(iV,:) = [XFine(Center(1)),YFine(Center(2))];
      case 'Average' % Compute a weighted average inside the fine beamforming 
        cP = cProjectionFine; 
        cP(cP<0) = 0; cP = cP.^2; cP = cP./sum(cP,'all');
        [cX,cY] = meshgrid(XFine,YFine);
        cR.Location(iV,:) = [sum(cP.*cX,'all'),sum(cP.*cY,'all')];
    end
    
    % ESTIMATE A MEASURE OF SNR / LOCATION CERTAINTY
    MAX = max(cProjectionRefined(:));
    MEDIAN = median(cProjectionRefined(:));
    STD = std(cProjectionRefined(:));
    cR.LocationCertainty(iV) = STD/(MAX-MEDIAN);
  
    % COMPUTE THE DIRECTIONALITY OF THE ELIPSE OF THE FINE PROJECTION
    % This assumes equal resolution in both directions
    if P.ComputeDirection
      ContourLevel = 0.8*(MAX-MIN)+MIN;
      C = contourc(cProjectionFine,[ContourLevel,ContourLevel]);
      NSteps = C(2,1);
      C = C(:,2:NSteps+1); % first row is X (Dim = 1 for Projection), second row is Y ()
      [DistMax,cIndMax] = max(sqrt(sum([C(1,:)-Center(1);C(2,:)-Center(2)].^2,1)));
      [DistMin,cIndMin] = min(sqrt(sum([C(1,:)-Center(1);C(2,:)-Center(2)].^2,1)));
      Vec = [C(1,cIndMax)-Center(1);C(2,cIndMax)-Center(2)]; % [X;Y]
      Vec = ((DistMax-DistMin)/DistMin-1)*Vec;
      cR.Vec(iV,1:2) = Vec;
    end

    % PLOT RESULTS
    if P.FIG
      figure(P.FIG); clf;[~,AH] =  axesDivide(2,1,'c');
      axes(AH(1)); hold on
      h=pcolor(Localization.X,Localization.Y,Localization.Data.Data.Frames(:,:,iV));
      set(h,'EdgeColor','None'); axis tight;
      axes(AH(2)); hold on
      h = pcolor(XIndFine,YIndFine,cProjectionFine);
      set(h,'EdgeColor','None'); axis tight;
      plot3([Center(1),Center(1)+Vec(1)],[Center(2),Center(2)+Vec(2)],[1,1],'k');
      plot3(C(1,:),C(2,:),repmat(MAX,[1,NSteps]));
      plot3(Center(1),Center(2),1,'r.')
      set(gca,'DataAspectRatio',[1,1,1]);
    end
  else
    cR.Vec(iV,1:2) = [NaN,NaN];
    cR.Location(iV,:) = [NaN,NaN];
    cR.LocationCertainty(iV) = NaN;
    continue;
  end
end
cR.Dims = size(cR.Frames);
Localization.Data.Data = cR;
