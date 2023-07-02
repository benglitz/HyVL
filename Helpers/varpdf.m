function Var = varpdf(x,px)

Var = px*((x-meanpdf(x,px)).^2)';