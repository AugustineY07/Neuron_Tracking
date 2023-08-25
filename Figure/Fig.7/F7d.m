
clear all

fig_path = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\Github Figure reproduce\F2'; %NEED CHANGE to local path
load(fullfile(fig_path,'F2d_count.mat'));

day = 5;
shank = 4;
dur = [1 2 13 23 48];


figure()
h = bar(KS_shankSum,0.5);
h.FaceColor = [0.30 0.75 0.93];
ylim([0 800])
xlabel('Days','FontSize',18,'FontWeight','Bold','Color','k')
ylabel('Number of Units across all shanks','FontSize',18,'FontWeight','Bold','Color','k')
title('AL032 Unit Counts','FontSize',18)



for id = 1:day
    xlb{id} = ['day ',num2str(dur(id))];
    text(id-0.2,KS_shankSum(id)/2,num2str(KS_shankSum(id)),'FontSize',13,'FontWeight','bold','HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'bottom') %KS text
end
for ir = 1:day-1
    if ir == 1
        r = max(KS_shankSum(1),KS_shankSum(ir+1));
    else
        r = max(KS_shankSum(1),KS_shankSum(2)) + (1000-700)/3*(ir-1);
    end
    b = drawbrace([1 r], [ir+1 r],12,'Color', 'k');
    text((ir+2)/2, 400+100*(ir-1),num2str(nRef_shankSum(ir)),'FontSize',13,'FontWeight','bold','HorizontalAlignment', 'center') %ref text
end
xticks([1:1:day])
xticklabels(xlb)
ax = gca;
ax.FontSize = 16; 
ax.Box = 'off';
set(ax,'TickDir','out');
l = [h,b];
legend(l,'KSgood units','reference pairs','Location','northeast')
