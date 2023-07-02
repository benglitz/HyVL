function [Data,F,T,Thresh,SRSpec] = HF_specgram(Data,Nfft,SRAI,FRange,NOverlap,Sparse,Abs)
% COMPUTE SPECTROGRAM WITH FREQUENCY RANGE REDUCTION

  switch NOverlap
    case Nfft/2
      SpecSteps = floor(length(Data)/Nfft*2);
      Data = reshape(Data(1:SpecSteps*Nfft/2),Nfft/2,SpecSteps);
      Data = [Data(:,1:end-1) ; Data(:,2:end)]; 
      SRSpec = SRAI/(Nfft/2);
    case 0
      SpecSteps = floor(length(Data)/Nfft);
      Data = reshape(Data(1:SpecSteps*Nfft),Nfft,SpecSteps);
      SRSpec = SRAI/Nfft;
    otherwise
      StepSize = Nfft-NOverlap;
      Steps = [1:StepSize:length(Data)-Nfft+1];
      SpecSteps = length(Steps);
      DataN = single(zeros(Nfft,SpecSteps));
      for i=1:SpecSteps
        DataN(:,i) = Data(StepSize*i:StepSize*i+Nfft-1);
      end
      Data = DataN; clear DataN;     
  end
  SpecSteps = SpecSteps - 1;
  
  % WINDOW DATA
  Window = hanning(Nfft);
  Data = bsxfun(@times,Data,Window);
  
  % Fourier transform
  Data = fft(Data,Nfft);
  Data = Data(1:end/2,:);
  if Abs Data = abs(Data); end
  
  % COMPUTE THE FREQUENCY AND TIME VECTORS
  F = linspace(0,SRAI/2,Nfft/2+1); F = F(1:end-1);
  if nargout > 2;   T = [0:1/SRSpec:(SpecSteps-1)/SRSpec]; end
  
  if nargin > 3 && ~isempty(FRange)
    % SELECT RANGE OF FREQUENCIES
    Ind = logical((F>=FRange(1)) .* (F<=FRange(2)));
    F = F(Ind);
  
    % Assume input signal is real, and take only the necessary part
    Data = Data(Ind,:);
  end
  
  if Sparse
    % Convert to sparse representation based on threshold
    ThreshData = Data([round(end/2),round(2*end/3)],end-min(10000,size(Data,2)-1):end-5000);
    Thresh = prctile(ThreshData(:),70);
    Ind = find(abs(Data) >= Thresh);
    [Ind1,Ind2] = ind2sub(size(Data),Ind);
    Data= sparse(Ind1,Ind2,double(Data(Ind)),size(Data,1),size(Data,2));
  end