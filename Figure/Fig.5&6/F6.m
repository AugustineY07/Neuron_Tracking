% Unit 9 from AL032 Shank 1
clear all
close all
clc

addpath(genpath('D:\Data\Pipeline\npy\npy-matlab'))

fig_path = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\Github Figure reproduce\F5&6'; %NEED CHANGE to local path
load(fullfile(fig_path,'F6_data.mat'))


for id = 1:day-1
    % plot waveforms
    plot_2_wfs('all',peakWf_all,chan_pos,L2_weight_full,clu_in_chain1,clu_in_chain2,day,dur,id);
end

% plot locations
h2 = figure();
h2.Name = sprintf('Chain %d',2);
h2.Units = 'centimeters';
set(h2,'Position',[6.1 6.6 7.6 22.5])
scatter(chan_pos(:,1),chan_pos(:,2),100,[0.9290 0.6940 0.1250],'square','filled'); hold on
plot_loc(loc_ref,day);
for ll = 1:day; unit{ll} = sprintf('Dataset %d unit',ll); end
xlim([min(chan_pos(:,1))-100 max(chan_pos(:,1))+68])
ylim([min(locs(:,2))-20, max(locs(:,2))+20])
legend(['Electrodes',unit],'Location','west')

% plot vfp and firing rate
plot_vfp('reference',rsp_chain,psth_chain,dur,day,rsp,psth)
plot_fr('reference',fr_all,fr_change_full,day);














%--------------------------------------helper functions--------------------------------------
function plot_fr(type,f_r,fr_change,day)
figure()
subplot(2,1,1)
plot(f_r(:),'LineWidth',2)
title(sprintf('%s chain %d firing rate',type,2))
xticks([1:1:day])
ylim([0 max(f_r(:))+1])
ylabel('Firing rate, spikes/s','FontSize',18,'FontWeight','Bold','Color','k')
ax = gca;
ax.FontSize = 16; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside

subplot(2,1,2)
plot(fr_change(:),'LineWidth',2)
xticks([1:1:day-1])
ylabel('Normalized firing change compare to prior dataset','FontSize',18,'FontWeight','Bold','Color','k')
title(sprintf('%s chain %d firing rate percent change',type,2),'FontSize',18,'FontWeight','Bold','Color','k')
ax = gca;
ax.FontSize = 16; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside
end


function plot_vfp(type,rsp_chain,psth_chain,dur,day,rsp,psth)
figure()
plot([1:1:day-1],rsp(1,:),'LineWidth',2); hold on
plot([1:1:day-1],psth(1,:),'LineWidth',2); hold on
sim = rsp(1,:)+psth(1,:);
plot([1:1:day-1],sim,'LineWidth',2); hold on
yline(1,'k--','LineWidth',2); hold on %ref threshold
ylabel('correlations or correlation sum','FontSize',18,'FontWeight','Bold','Color','k')
legend('vfp','PSTH','simScore','simScore threshold')
title(sprintf('Subject %s Shank %d %s chain %d fingerprint','AL032',3,type,2),'FontSize',18,'FontWeight','Bold','Color','k')
ax = gca;
ax.FontSize = 16; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside

% plot lowest score receptive field, maybe also plot the non-ref day for mixed and putative
figure()
for id = 1:day-1
    Dblue = [0 0.45 0.74];
    Dred = [0.64 0.08 0.18];

    rsp1 = rsp_chain{id};
    rsp2 = rsp_chain{id+1};
    psth1 = psth_chain{id};
    psth2 = psth_chain{id+1};

    subplot(2,day-1,id)
    [hAx,hLine1,hLine2] = plotyy(1:numel(rsp1),rsp1, 1:numel(rsp2),rsp2+max(rsp1));
    set(hAx,'FontSize', 18); %tick font
    set(hAx,'Box', 'off'); %remove tick box
    set(hAx,'TickDir','out'); %tickmark towards outside
    p1 = get(hAx,'position');
    hLine1.LineWidth = 2;
    hLine1.Color = Dblue;
    hLine2.LineWidth = 2;
    hLine2.Color = Dred;
    set(hAx,'ylim',[0 max(rsp2)+max(rsp1)+5])
    set(hAx(1),'YTick',[0:5:max(rsp2)+max(rsp1)])
    set(hAx(1),'YColor',Dblue);
    set(hAx(2),'YTick',[max(rsp1):5:max(rsp1)+15])
    set(hAx(2),'yticklabels',[0:5:max(rsp2)+max(rsp1)+15])
    set(hAx(2),'YColor',Dred);
    xlabel('Stimulus id','FontSize',18,'FontWeight','Bold','Color','k')
    ylabel(sprintf('Spike rate (day %d)',dur(id)),'FontSize',18,'FontWeight','Bold','Color',Dblue) % left y-axis
    ylabel(hAx(2),sprintf('Spike rate (day %d)',dur(id+1)),'FontSize',18,'FontWeight','Bold','Color',Dred) % right y-axis
    legend({sprintf('day %d',dur(id)),sprintf('day %d',dur(id+1))},'Location','northeast','Orientation','vertical','FontSize',18)
    h1 = subtitle(sprintf('corr\\_vfp = %.2f', rsp(id)));
    h1.FontSize = 18;

    subplot(2,day-1,id+day-1)
    [hAx2,hLine3,hLine4] = plotyy(1:numel(psth1),psth1, 1:numel(psth2),psth2+max(psth1));
    set(hAx2,'FontSize', 18); %tick font
    set(hAx2,'Box', 'off'); %remove tick box
    set(hAx2,'TickDir','out'); %tickmark towards outside
    hLine3.LineWidth = 2;
    hLine3.Color = Dblue;
    hLine4.LineWidth = 2;
    hLine4.Color = Dred ;
    set(hAx2,'ylim',[0 max(psth2)+max(psth1)+0.4])
    set(hAx2(1),'YTick',[0:0.2:max(psth2)+max(psth1)])
    set(hAx2,'XTick',[0:500:2000])
    set(hAx2(1),'YColor',Dblue);
    set(hAx2(2),'YTick',[max(psth1):0.2:max(psth1)+1])
    set(hAx2(2),'yticklabels',[0:0.2:max(psth2)+max(psth1)+1])
    set(hAx2(2),'YColor',Dred);
    xlabel('Time (ms)','FontSize',18,'FontWeight','Bold','Color','k')
    ylabel(sprintf('Normalized spike rate (day %d)',dur(id)),'FontSize',18,'FontWeight','Bold','Color',Dblue) % left y-axis
    ylabel(hAx2(2),sprintf('Normalized spike rate (day %d)',dur(id+1)),'FontSize',18,'FontWeight','Bold','Color',Dred) % right y-axis
    legend({sprintf('day %d',dur(id)),sprintf('day %d',dur(id+1))},'Location','northeast','Orientation','vertical','FontSize',18)
    h2 = subtitle(sprintf('corr\\_PSTH = %.2f', psth(id)));
    h2.FontSize = 18;
    ax = gca;
    ax.FontSize = 16; %tick font
    ax.Box = 'off'; %remove tick box
    set(ax,'TickDir','out'); %tickmark towards outside
end

end


function plot_loc(locations,day)
for id = 1:day-1
    cmap = colormap(hsv);
    c = cmap(4+42*((1:day)-1),:);
    scatter(locations{1,2}(id,1),locations{1,2}(id,2),50,c(id,:),'o','filled'); hold on
    if id == day-1
        scatter(locations{1,2}(id+1,1),locations{1,2}(id+1,2),50,c(id+1,:),'o','filled'); hold on
    end
end
ax = gca;
ax.FontSize = 16; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside
end




