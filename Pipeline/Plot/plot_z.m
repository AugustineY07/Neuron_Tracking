function plot_z(path)
% plot z distance distribution of all matched units

load(fullfile(path,'Output.mat'));
matched_units = output.results_wth;
all_z = matched_units(:,7);

figure()
histogram(all_z,20)
xlabel('z distance (\mum)')
ylabel('Number of units')