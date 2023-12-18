function pdelta = pd_fg_est_user(delta, a,b,c,d)

% Fit the distribution of all pairs, with no z distance threshold, from an
% emd matching algorithm

% In matlab, can express as a surface fit (rather than minimization)
% 'Fit objects' in MATLAB order coeffcients alphabetically, so they are
% nicknamed (a,b,c) for the input

% distribution is estiamted as the sum of two exponentials, for correct and
% incorrect pairs. Decay constant (tau) assumed known for the correct pairs

corr_sigma = a;  % width of distribution of true pairs, known from reference pairs
normFactor = b;
incorr_tau = 1/c;
frac_corr = d;

fg_norm = 2/(corr_sigma*sqrt(2*pi));

pdelta = normFactor*( frac_corr*fg_norm*exp(-(delta.^2/(2*corr_sigma^2))) + ...
                             (1-frac_corr)*incorr_tau*exp(-incorr_tau.*delta) );



end