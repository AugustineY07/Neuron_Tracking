clear all

fig_path = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\Github Figure reproduce\F3&4'; %NEED CHANGE to local path
load(fullfile(fig_path,'F3_data.mat'))

% set z_dist value threshold array
dist = EMD_pair(:,7);
dist_ref = all_pair_ref(:,7);
thre_max = 100;
thre_min = min(dist_ref);
threshold = [ceil(thre_min):1:thre_max];
threshold_summary = [];
for ith = 1:length(threshold)
    % loop: find ref_pairs in selected subset
    th = threshold(ith);
    idx_subth = find(dist<=th);
    KSgood_in_range = EMD_pair(idx_subth,:);
    ref_in_range = KSgood_in_range(KSgood_in_range(:,1)==1,:);
    num_ref_correct = sum(ref_in_range(:,8)); %num correct pairs in in-threshold pairs
    num_ref_all = size(ref_in_range,1);
    num_KS = size(KSgood_in_range,1);
    accuracy = num_ref_correct/num_ref_all;
    threshold_summary(ith,:) = [th; num_ref_correct; num_KS; accuracy];
end %end of threshold loop

% find 95% threshold
th_95 = find(threshold_summary(:,4) >= 0.95, 1, 'last');

% plot distance_threshold v.s. accuracy for all data
figure()
[ax, p1, p2]=plotyy(threshold,threshold_summary(:,4),threshold,threshold_summary(:,3));
p1.Color = [0 0.45 0.74];
p1.LineWidth = 3;
p2.Color = [0.64 0.08 0.30];
p2.LineWidth = 3;
hold(ax(1),'on')
h1 = scatter(ax(1),threshold,threshold_summary(:,4),'.','MarkerFaceColor',[0 0.45 0.74]);
hold(ax(2),'on')
h2 = scatter(ax(2),threshold,threshold_summary(:,3),'.','MarkerFaceColor',[0.64 0.08 0.30]);
ref_idx = find(EMD_pair(:,1) == 1);
z_mean = mean(EMD_pair(ref_idx,7)); %mean z-dist
h4 = xline(repelem(z_mean,2),'k','LineWidth',1);
h5 = xline(10,'k--','LineWidth',1);
set(ax(1),'YLim',[0.8 1])
set(ax(2),'YLim',[100 1010])
set(ax(1),'YTick',[0.84:0.02:1],'FontSize',24,'YColor',[0 0.45 0.74])
set(ax(2),'YTick',[100:100:1000],'FontSize',24,'YColor',[0.64 0.08 0.30])
set(ax(2),'TickDir','out');
xlabel('Threshold z-distance','FontSize',30,'FontWeight','Bold','Color','k')
ylabel(ax(1),'Accuracy of ref pairs included','FontSize',30,'FontWeight','Bold','Color',[0 0.45 0.74])
ylabel(ax(2),'Number of ref pairs included','FontSize',30,'FontWeight','Bold','Color',[0.64 0.08 0.30])
h = [p1(1);p2;h4;h5]; 

legend(h,{'sampled accuracy','KS pair count','mean ref Z-dist','','threshold at Z = 10\mum'},'Location','best','FontSize',30); 
legend boxoff 
title('All sampled accuracy v.s. threshold','FontSize',30)
ax = gca;
ax.FontSize = 24; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out');



