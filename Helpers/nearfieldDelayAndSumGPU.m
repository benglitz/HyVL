function SteeredPower = nearfieldDelayAndSumGPU(Data,Xm,Ym,Zm,Xb,Yb,Distance,SpeedOfSound,SR)
% nearfieldDelayAndSum - calculate delay and sum in time domain

Data = gpuArray(single(Data));
Xm = gpuArray(single(Xm));
Ym = gpuArray(single(Ym));
Zm = gpuArray(single(Zm));
Xb = gpuArray(single(Xb));
Yb = gpuArray(single(Yb));
Distance = gpuArray(single(Distance));
Scaler = gpuArray(single(SpeedOfSound/SR));

NSteps = size(Data,1);
NMic = size(Data,2);
SteeredPower = zeros(length(Xb),length(Yb),'gpuArray');
for iX = 1:length(Xb) % Loop over Scanpoints in X
  for iY = 1:length(Yb) % Loop over Scanpoints in Y
    
    % Distance from scanning point to microphones
    cXb = Xb(iX); cYb = Yb(iY);
    Distances = sqrt((Xm - cXb).^2 + (Ym - cYb).^2 + (Zm - Distance).^2);
    Distances = Distances - min(Distances); % Makes Distances Start from 0, but keeps relative relation
    
    %Delay and sum
    dt = Distances/Scaler; %time delay in samples, not rounded to perform subsampling
    dtMax = ceil(max(dt));
    dtR = round(dt);
    cNSteps = NSteps-dtMax-1;
    cInd = [1:cNSteps]';
    Inds = (dtR + cInd) + ((NSteps*ones(cNSteps,1,'gpuArray'))*[0:NMic-1]);
    SteeredResp = mean(Data(Inds), 2); % Sum across microphones (keeping NANs, to avoid summing over subset of microphones)
    SteeredPower(iX,iY) = sum(SteeredResp.^2)/(cNSteps); % Now only summing across the time samples that are not NaN 
  end
end
SteeredPower = gather(SteeredPower);