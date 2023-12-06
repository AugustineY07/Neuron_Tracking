function plot_unit(chan_pos,input,result_file,EMD_data)
% plot matched units on a probe

load(fullfile(input.input_path,result_file,'Output.mat'));
matched_units = output.results_wth;
points = load(fullfile(input.EMD_path,EMD_data));
f1 = points.f1;
f1_label = points.f1_labels;
f2 = points.f2;
f2_label = points.f2_labels;



h1 = figure();
% h1.Position = [1680,41,305,1140];
scatter(chan_pos(:,1),chan_pos(:,2),100,[0.9290 0.6940 0.1250],'square','filled'); hold on
for ip = 1:size(matched_units,1)
    unit1 = matched_units(ip,3);
    unit2 = matched_units(ip,2);
    idx_label1 = find(f1_label==unit1);
    idx_label2 = find(f2_label==unit2);
    loc_unit1 = f1(idx_label1,1:3);
    loc_unit2 = f2(idx_label2,1:3);
    c = rand(1, 3);
    scatter(loc_unit1(:,1), loc_unit1(:,2),50,c,'o'); hold on;
    scatter(loc_unit2(:,1), loc_unit2(:,2),50,c,'o','filled'); hold on;
    u = loc_unit2(:,1)-loc_unit1(:,1);
    v = loc_unit2(:,2)-loc_unit1(:,2);
    quiver(loc_unit1(:,1),loc_unit1(:,2),u,v,'Color','black','LineStyle',':','LineWidth', 1,'MaxHeadSize', 1);
end
xlim([min(chan_pos(:,1))-20,max(chan_pos(:,1))+20]);
ylim([min(chan_pos(:,2))-50,max(chan_pos(:,2))+50]);
% ylim([4500,5000]);
axis([min(chan_pos(:,1))-20, max(chan_pos(:,1))+20, [min(chan_pos(:,2))-50,max(chan_pos(:,2))+50]]);
axis equal




