logit(x::Float64) = log(x/(1-x))
invlogit(x::Float64) = exp(x)/(1+exp(x))
