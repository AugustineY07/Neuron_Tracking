function fit_zDist

clear all

addpath(genpath('C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\')) % NEED CHANGE: Path to you code

% load in arrays of measured distances
dataDir = 'Z:\zDist_distribution_fitting\'; % NEED CHANGE: Path to you data
load(fullfile(dataDir,"KS_all.mat")) % NEED CHANGE: Data file name
zDist_all = EMD_pair(:,7)'; % NEED CHANGE: z distances

% edges for histogramming the distributions
edges = (0:2:100);
bin_centers = (edges(1:numel(edges)-1) + 1)';

all_cnts = histcounts(zDist_all, edges);




% fit without reference
fitP_FG = fittype('pd_fg_est_user(delta, a, b, c, d)', ...
    'independent', 'delta', 'dependent', 'pdelta');

a0 = 3.3837;
b0 = sum(all_cnts);  % est norm factor from integral of all_cnts
c0 = 50; % z_half, estimate
d0 = 0.5; % fraction correct, estimate

up_bound = [inf, inf, inf, 1];
low_bound = [-inf, -inf, -inf, 0];

all_fo_start = fit( bin_centers, all_cnts', fitP_FG, 'StartPoint', [a0, b0, c0, d0], ...
              'upper', up_bound, 'lower', low_bound)

CIF = predint(all_fo_start,bin_centers,0.95,'Functional');


% Plot the original data and the fitted curve
figure();
% plot((1:length(all_cnts))', all_cnts, '-', 'DisplayName', 'all_cnts');
pdelta_calc = feval(all_fo_start,bin_centers);
plot(bin_centers,pdelta_calc,'LineWidth',3); hold on;
scatter( bin_centers, all_cnts,30, 'filled', 'square'); hold on;
plot(bin_centers,CIF,':k','LineWidth',1.5)
ylim([0 inf])
title('Fit Result for z Distance distributions');
legend('Location', 'Best');
legend('Folded Gaussian + Exponential fit','Unit count data','95% CI upper bound','95% CI lower bound')
xlabel('Z distances (\mu m)');
ylabel('Counts');
ax = gca; 
ax.FontSize = 18;
hold off;
ax.Box = 'off';
set(ax,'TickDir','out');






% Plot fraction fp vs threshold 
figure()
th_vals = (0:2:20);

decay = all_fo_start.a;
norm_factor = all_fo_start.b;
frac_corr_fit = all_fo_start.d;
incorr_zhalf = all_fo_start.c;
incorr_tau = 1/incorr_zhalf;
hits_th = norm_factor * frac_corr_fit * normcdf(th_vals,0,frac_corr_fit);%frac_corr_fit
fp_th = norm_factor * (1-frac_corr_fit) * (1-exp(-incorr_tau.*th_vals));

xRmin = 10;
yRmin = (norm_factor * (1-frac_corr_fit) * (1-exp(-incorr_tau.*xRmin))) ./ (norm_factor * frac_corr_fit * normcdf(xRmin,0,frac_corr_fit) + norm_factor * (1-frac_corr_fit) * (1-exp(-incorr_tau.*xRmin)));
plot(th_vals, fp_th./(hits_th + fp_th),'LineWidth',3); hold on;
plot([0 xRmin], [1 1]*yRmin, '--r'); hold on;
plot([1 1]*xRmin, [0 1]*yRmin, '--r'); hold on;
plot(xRmin, yRmin, 'or'); hold on;
text(0,yRmin, sprintf('%.2f',yRmin), 'Horiz','right', 'Vert','middle','FontSize',16,'Color','r')     % Horizontal Line Axis Label

ylim([0 inf])
title('False Positive fraction vs. \Deltaz threshold','FontSize',18,'FontWeight','Bold','Color','k');
ylabel('Fraction false positive','FontSize',18,'FontWeight','Bold','Color','k');
xlabel('\Deltaz threshold (\mum)','FontSize',18,'FontWeight','Bold','Color','k');
ax = gca;
ax.FontSize = 18; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside

end

