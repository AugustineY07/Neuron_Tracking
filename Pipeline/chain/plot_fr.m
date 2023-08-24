function fr_change = plot_fr(fr_all, fr_change, numData, ichain)

figure()
subplot(2,1,1)
plot(fr_all(ichain,:),'LineWidth',2)
title(sprintf('Chain %d firing rate',ichain))
xticks([1:1:numData])
ylim([0 max(fr_all(ichain,:))+1])
ylabel('Firing rate, spikes/s','FontSize',18,'FontWeight','Bold','Color','k')
ax = gca;
ax.FontSize = 16; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside

subplot(2,1,2)
plot(fr_change(ichain,:),'LineWidth',2)
xticks([1:1:numData-1])
ylabel('Normalized firing change compare to prior dataset','FontSize',18,'FontWeight','Bold','Color','k')
title(sprintf('Chain %d firing rate percent change',ichain),'FontSize',18,'FontWeight','Bold','Color','k')
ax = gca;
ax.FontSize = 16; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside
end





