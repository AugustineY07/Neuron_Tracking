function pdelta = FP_exp_user(delta, a, b)

decay = a;
normFactor = b;

pdelta = normFactor*exp(-decay*delta);
end