function plot_z(path)
% plot z distance distribution of all matched units

load(fullfile(path,'Output.mat'));
matched_units = output.results_wth;
all_z = matched_units(:,7);
unthreshold_matched_units = output.all_results_post;
unthreshold_z = unthreshold_matched_units(:,7);

figure()
subplot(2,1,1)
histogram(all_z,20)
xlim([0 10])
xlabel('z distance (\mum)')
ylabel('Number of units')
subtitle('threshold z distribution')

subplot(2,1,2)
histogram(unthreshold_z,70)
xlim([0 100])
xlabel('z distance (\mum)')
ylabel('Number of units')
subtitle('all z distribution')