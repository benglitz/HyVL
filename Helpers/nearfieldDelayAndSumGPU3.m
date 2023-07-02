function SteeredPower = nearfieldDelayAndSumGPU3(Data,Xm,Ym,Zm,Xb,Yb,Distance,SpeedOfSound,SR)
% nearfieldDelayAndSum - calculate delay and sum in time domain

Data = gpuArray(single(Data));
Xm = gpuArray(single(Xm));
Ym = gpuArray(single(Ym));
Zm = gpuArray(single(Zm));
Xb = gpuArray(single(Xb));
Yb = gpuArray(single(Yb));
Distance = gpuArray(single(Distance));
Scaler = gpuArray(single(SpeedOfSound/SR));

NSteps = gpuArray(size(Data,1));
NMic = gpuArray(size(Data,2));
SteeredPower = zeros(length(Xb),length(Yb),'gpuArray');
Tmp = NSteps*ones(NSteps,1,'gpuArray')*[0:NMic-1];
[Xb,Yb] = meshgrid(Xb,Yb);
Xb = permute(Xb,[3,1,2]);
Yb = permute(Yb,[3,1,2]);

% Distance from scanning point to microphones
Distances = sqrt((Xm' - Xb).^2 + (Ym' - Yb).^2 + (Zm' - Distance).^2);
Distances = Distances - min(Distances); % Makes Distances Start from 0, but keeps relative relation

%Delay and sum
dt = Distances/Scaler; %time delay in samples, not rounded to perform subsampling
dtMax = ceil(max(dt(:)));
dtR = round(dt);
cNSteps = NSteps-dtMax-1;
cInd = [1:cNSteps]';
Inds = (permute(dtR,[4,1,2,3]) + cInd) + Tmp(1:cNSteps,:);
SteeredResp = mean(Data(Inds), 2); % Sum across microphones
SteeredPower = sum(SteeredResp.^2,1)./cNSteps; % Now only summing across the time samples that are not NaN
SteeredPower = gather(SteeredPower);