function estFP_10um = fit_zDist( varargin )

if (length(varargin) == 0)
    clear all;
    % Edit parameters manually
    % Input data for the fit is an array of z distances, here extracted
    % from a summary datafile. In standard pipeline output, z distances are
    % in column 7 of EMD_postN.all_results or in the results folder Output.mat 
    % output.all_results_post (along with more diagnostic metrics).
    genpath('C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\') % NEED CHANGE: Path to the repo
    % load in arrays of measured distances. This can be any array in a mat
    % file, replace these lines to load in your data. Note that the z
    % distances to are recorded in column all_results(:,7) of the
    % EMD_postN.mat file output for each match. To determine a good value
    % for the width of the folded gaussian distribution of correct matches,
    % pool results from several low drift matched datasets.
    zDistMat = load('zDist_all.mat');
    zDist_all = zDistMat.data; % NEED CHANGE: z distances
    data_label = '';
    % set fg_width -1 to fit the width of the folded gaussian 
    % of correct pairs; set to > 0 to fix at a known width; 3.5 typical for
    % NP 2.0 data
    fg_width = -1; 
else
    % fit the zDist_all data passed in 
    inputCell = varargin(1);
    zDist_all = inputCell{1};   % array of distances
    inputCell = varargin(2);
    data_label = inputCell{1};  % label for graphs
    % 3rd argument: set to -1 to fit the width of the folded gaussian 
    % of correct pairs; set to > 0 to fix at a known width; 3.5 typical for
    % NP 2.0 data
    inputCell = varargin(3);
    fg_width = inputCell{1};  % as above, call with -1 to fit fg_width, or a value > 0 to fix.
end 

% edges for histogramming the distributions
edges = (0:2:100);
bin_centers = (edges(1:numel(edges)-1) + 1)';

all_cnts = histcounts(zDist_all, edges);



fitP_FG = fittype('pd_fg_est_user(delta, a, b, c, d)', ...
        'independent', 'delta', 'dependent', 'pdelta');

if fg_width < 0
    % fit with no assumption for the width of the folded Gaussian distribtuion
    % of correct pairs
    a0 = 3.3837;
    b0 = sum(all_cnts);  % est norm factor from integral of all_cnts
    c0 = 50; % z_half, estimate
    d0 = 0.5; % fraction correct, estimate
    
    up_bound = [inf, inf, inf, 1];
    low_bound = [-inf, -inf, -inf, 0];
    
    all_fo_start = fit( bin_centers, all_cnts', fitP_FG, 'StartPoint', [a0, b0, c0, d0], ...
                  'upper', up_bound, 'lower', low_bound)
    
    CIF = predint(all_fo_start,bin_centers,0.95,'Functional');
else
    % fix width of correct pairs
    
    a0 = fg_width;
    b0 = sum(all_cnts);  % est norm factor from integral of all_cnts
    c0 = 50; % z_half, estimate
    d0 = 0.5; % fraction correct, estimate
    
    up_bound = [inf, inf, 1];
    low_bound = [10, 0, 0];

    fitP_FG_fix = fittype('pd_fg_est(delta, a, b, c, d)', ...
    'independent', {'delta'}, 'dependent', 'pdelta', 'problem','a' );

    all_fo_start = fit( bin_centers, all_cnts', fitP_FG_fix, 'problem', a0, 'StartPoint', [b0, c0, d0], ...
              'upper', up_bound, 'lower', low_bound)
  
    CIF = predint(all_fo_start,bin_centers,0.95,'Functional');
end


% Plot the original data and the fitted curve
figure();
% plot((1:length(all_cnts))', all_cnts, '-', 'DisplayName', 'all_cnts');
pdelta_calc = feval(all_fo_start,bin_centers);
plot(bin_centers,pdelta_calc,'LineWidth',3); hold on;
scatter( bin_centers, all_cnts,30, 'filled', 'square'); hold on;
plot(bin_centers,CIF,':k','LineWidth',1.5)
ylim([0 inf])
if isempty(data_label)
    title('Fit Result for z Distance distributions');
else
    title(sprintf('Fit Result for %s', data_label),'Interpreter','none');
end
legend('Location', 'northeast');
legend('Folded Gaussian + Exponential fit','Unit count data','95% CI upper bound','95% CI lower bound')
xlabel('Z distances (\mu m)');
ylabel('Counts');
ax = gca; 
ax.FontSize = 12;
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
if isempty(data_label)
    title('False Positive fraction vs. \Deltaz threshold','FontSize',18,'FontWeight','Bold','Color','k');
else
    title(sprintf('FP vs. thresh for %s', data_label),'Interpreter','none');
end
ylabel('Fraction false positive','FontSize',18,'FontWeight','Bold','Color','k');
xlabel('\Deltaz threshold (\mum)','FontSize',18,'FontWeight','Bold','Color','k');
ax = gca;
ax.FontSize = 18; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside


estFP_10um = fp_th(th_vals==10);

end

