function fit_zDist_exp

close all
clear all

fig_path = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\Github Figure reproduce\F3&4'; %NEED CHANGE to local path
load(fullfile(fig_path,'F4a_data.mat'))



% edges for histogramming the distributions
edges = (0:2:100);
bin_centers = (edges(1:numel(edges)-1) + 1)';

% fit hits with folded gaussian
fitFG = fittype('folded_gauss(delta, a, b)', ...
    'independent', {'delta'}, 'dependent', 'pdelta' );

% fit all units
fitP_FG = fittype('pd_fg_est(delta, a, b, c, d)', ...
    'independent', {'delta'}, 'dependent', 'pdelta', 'problem','a' );

% fit FP
fitFP = fittype('FP_exp(delta, a, b)', ...
    'independent', {'delta'}, 'dependent', 'pdelta' );



ref_cnts = histcounts(zDist_ref, edges);
all_cnts = histcounts(zDist_all, edges);
fp_cnts = histcounts(zDist_incorr, edges);

% fit the distribution of ref pairs to a folded gaussian

a0 = 5;  %estimated, from falloff of the histogram
b0 = sum(ref_cnts); %normalization factor is the integral over all bins
up_bound = [2*a0, 2*b0];
low_bound = [0.2*a0, 0.2*b0];

ref_fo = fit( bin_centers, ref_cnts', fitFG, 'StartPoint', [a0, b0], ...
              'upper', up_bound, 'lower', low_bound)



% a2. plot fit to reference pairs
figure()
pdelta = feval(ref_fo,bin_centers);
plot(bin_centers,pdelta,'LineWidth',3); hold on %plot(ref_fo,bin_centers,ref_cnts)
scatter( bin_centers, ref_cnts,30, 'filled', 'square');hold on;
title('Ref pairs','FontSize',18,'FontWeight','Bold','Color','k');
ylabel('Number of units','FontSize',18,'FontWeight','Bold','Color','k');
xlabel('\Deltaz (\mum)','FontSize',18,'FontWeight','Bold','Color','k');
legend('Folded Gaussian fit','Unit count data')
ax = gca;
ax.FontSize = 16; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside





% fit distribution of all pairs

ref_sigma = ref_fo.a

b0 = sum(all_cnts);  % est norm factor from integral of all_cnts
c0 = 70; % z_half, estimate
d0 = 0.5; % fraction correct, estimate

up_bound = [b0*3, c0*3, 0.9];
low_bound = [b0/3, c0/3, 0.1];

all_fo = fit( bin_centers, all_cnts', fitP_FG, 'problem', ref_sigma, 'StartPoint', [b0, c0, d0], ...
              'upper', up_bound, 'lower', low_bound)
    


% a1. Plot sums as points + fits 
figure();
pdelta_calc = feval(all_fo,bin_centers);
plot(bin_centers,pdelta_calc,'LineWidth',3); hold on;
scatter( bin_centers, all_cnts,30, 'filled', 'square'); hold on;
title('All pairs','FontSize',18,'FontWeight','Bold','Color','k');
ylabel('Number of units','FontSize',18,'FontWeight','Bold','Color','k');
xlabel('\Deltaz (\mum)','FontSize',18,'FontWeight','Bold','Color','k');
legend('Folded Gaussian+Exponential fit','Unit count data')
ax = gca;
ax.FontSize = 16; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside




% a3. Plot fraction fp vs threshold 
figure()
th_vals = (0:2:20);

decay = ref_fo.a;
norm_factor = all_fo.b;
frac_corr_fit = all_fo.d;
incorr_zhalf = all_fo.c;
incorr_tau = 1/incorr_zhalf;
hits_th = norm_factor * frac_corr_fit * normcdf(th_vals,0,ref_sigma);
fp_th = norm_factor * (1-frac_corr_fit) * (1-exp(-incorr_tau.*th_vals));

plot(th_vals, fp_th./(hits_th + fp_th),'LineWidth',3);

title('False Positive fraction vs. \Deltaz threshold','FontSize',18,'FontWeight','Bold','Color','k');
ylabel('Fraction false positive','FontSize',18,'FontWeight','Bold','Color','k');
xlabel('\Deltaz threshold (\mum)','FontSize',18,'FontWeight','Bold','Color','k');
ax = gca;
ax.FontSize = 16; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside


end