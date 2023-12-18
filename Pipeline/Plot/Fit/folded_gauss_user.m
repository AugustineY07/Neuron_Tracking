function pdelta = folded_gauss_user(delta, a, b)

% For a range of delta >=0 , calculate the folded gaussian distribution
% with mean = 0, sigma = a, normalization = b
% 'Fit objects' in MATLAB order coeffcients alphabetically, so they are
% nicknamed (a,b,c) for the input


sigma = a;
normFactor = b;

fg_norm = 2/(sigma*sqrt(2*pi));

pdelta = normFactor*fg_norm*exp(-(delta.^2/(2*sigma^2)));

end