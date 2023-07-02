function SteeredPower = nearfieldDelayAndSum(Data,Xm,Ym,Zm,Xb,Yb,Distance,SpeedOfSound,SR)
% nearfieldDelayAndSum - calculate delay and sum in time domain

[NSteps,NMic] = size(Data);
assert(NMic<=64,'Check Data Dimensions!');
NXb = length(Xb);
NYb = length(Yb);

Data = double(Data);
SteeredPower = zeros(NXb,NYb);
Offsets = repmat([0:NSteps:(NMic-1)*NSteps],NSteps,1);
TBase = [0:NSteps-1];

parfor iX = 1:NXb % Loop over Scanpoints in X
  for iY = 1:NYb % Loop over Scanpoints in Y
    
    % Distance from scanning point to microphones
    cXb = Xb(iX); cYb = Yb(iY);
    Distances = sqrt((Xm - cXb).^2 + (Ym - cYb).^2 + (Zm - Distance).^2);
    Distances = Distances - min(Distances); % Makes Distances Start from 0, but keeps relative relation
    
    %Delay and sum
    dt = Distances/SpeedOfSound*SR; %time delay in samples, not rounded to perform subsampling
    dtMax = ceil(max(dt));
    dtR = round(dt);
    cInd = [1:NSteps-dtMax]';
    Data_aligned = NaN*zeros(NSteps,NMic);
    for iM = 1:NMic
      %Data_aligned(1:NSteps-dtMax:iM) = Data(dtR(iM)+cInd,iM); % About 20% faster, results very similar
      Data_aligned(:,iM) = interp1(TBase,Data(:,iM),TBase + dt(iM),'pchip',NaN);
    end
    SteeredResp = mean(Data_aligned, 2); % Sum across microphones (keeping NANs, to avoid summing over subset of microphones)
    SteeredPower(iX,iY) = nansum(SteeredResp.^2)/sum(~isnan(SteeredResp)); % Now only summing across the time samples that are not NaN
  end
end