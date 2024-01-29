function plot_dist(input)
% plot waveform vs physical distance 

load(fullfile(input.result_path, 'Output.mat'));
matched_units = output.results_wth;
dist = matched_units(:,5:6);

figure()
scatter(dist(:,1),dist(:,2),10,'filled')
xlabel('Physical distance (\mum)')
ylabel('Waveform distance')
