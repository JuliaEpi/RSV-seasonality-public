gompertz(a, x) = a * exp(log(1/a) * (1 - exp(-x)))
scaled_logistic(a, x) = a * StatsFuns.logistic(x)
scaled_exp(a, x) = a * exp(x)
