function S = skewnesspdf(x,px)

S = px*((x-meanpdf(x,px)).^3)'/stdpdf(x,px)^3;