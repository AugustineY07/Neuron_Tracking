function plot_loc(all_input, x_loc_all,z_loc_all, chan_pos, numData,ichain)

h2 = figure();
h2.Units = 'centimeters';
set(h2,'Position',[6.1 6.6 7.6 22.5])
scatter(chan_pos(:,1),chan_pos(:,2),100,[0.9290 0.6940 0.1250],'square','filled'); hold on
for ll = 1:numData; unit{ll} = sprintf('Dataset %d unit',ll); end
xlim([min(chan_pos(:,1))-100 max(chan_pos(:,1))+68])
ylim([min(min(z_loc_all))-20, max(max(z_loc_all))+20])
legend(['Electrodes',unit],'Location','west')

cmap = colormap(hsv);
for id = 1:numData-1
    c = cmap(4+42*((1:numData)-1),:);
    scatter(x_loc_all(ichain,id),z_loc_all(ichain,id),50,c(id,:),'o','filled'); hold on
    if id == numData-1
        scatter(x_loc_all(ichain,id+1),z_loc_all(ichain,id+1),50,c(id+1,:),'o','filled'); hold on
    end
end
ax = gca;
ax.FontSize = 16; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside


end
