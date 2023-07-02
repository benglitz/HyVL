function DataCorrected = C_correctImageDistortion(Data,Strength,Zoom,Mode,DataCenter)

if nargin<4; Mode = 'Image'; end

switch Mode
  case 'Image'
    Image = Data;
    HalfY = size(Image,1)/ 2;
    HalfX = size(Image,2)/ 2;
    
    CorrectionRadius = sqrt(HalfY ^ 2 + HalfX ^ 2) / Strength; 
    
    DataCorrected = zeros(size(Image));
    
    for iX = 1:size(Image,2)
      for iY = 1:size(Image,1)
        NewX = iX - HalfX;
        NewY = iY - HalfY;
        
        Dist = sqrt(NewX ^ 2 + NewY ^ 2);
        R = Dist / CorrectionRadius;
        
        if R == 0
          Theta = 1;
        else
          Theta = atan(R) / R;
        end
        
        SourceX = round(HalfX + Theta * NewX * Zoom);
        SourceY = round(HalfY + Theta * NewY * Zoom);
        
        DataCorrected(iY,iX) = Image(SourceY,SourceX);
      end
    end
    
  case 'Markers'
    % Data needs to be relative to the center of the image
    % Eq.4 from https://www.spiedigitallibrary.org/journals/Optical-Engineering/volume-56/issue-01/013108/Correction-of-image-radial-distortion-based-on-division-model/10.1117/1.OE.56.1.013108.full?SSO=1
    CenterDist = sqrt(sum((Data - DataCenter).^2,2));
    Lambda = Strength;
    DataCorrected = Data.*(1+Lambda*CenterDist.^2);
    
end