function S = kurtosispdf(x,px)

S = px*((x-meanpdf(x,px)).^4)'/stdpdf(x,px)^4;