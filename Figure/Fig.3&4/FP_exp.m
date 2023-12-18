function pdelta = FP_exp(delta, a, b)

decay = a;
normFactor = b;

pdelta = normFactor*exp(-decay*delta);
end