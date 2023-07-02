function Bins = HF_autobin(Min,Max,N)

Step = (Max - Min)/(N-1);
DecPos = log10(Step);
cDecPos = ceil(-DecPos);
cStep = 10^cDecPos*Step;
cStep = round(cStep);
Step = cStep * 10^-cDecPos;

cMin = Step*floor(Min/Step);
cMax = Step*ceil(Max/Step);

Bins = [cMin:Step:cMax];