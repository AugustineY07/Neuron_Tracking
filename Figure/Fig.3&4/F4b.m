% plot stacked bar accuracy results

clear all

fig_path = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\Github Figure reproduce\F3&4'; %NEED CHANGE to local path
data = load(fullfile(fig_path,'F4b_data.mat'));
fn = fieldnames(data);

subject = {'AL031','AL032','AL036'};%,'AL032','AL036'


for is = 1:length(subject)
    switch subject{is}
        case 'AL036'
            day = 6;
            shank = 4;
            dur = [1 2 4 10 16 24 36];
        case 'AL031'
            day = 4;
            shank = 1;
            dur = [1 2 11 36 45];
        case 'AL032'
            day = 4;
            shank = 4;
            dur = [1 2 13 23 48];
    end

    for id = 1:day
        xlb{id} = ['day 1-',num2str(dur(id+1))];
        groupLabels{id} = xlb{id};
    end





             %correct: blue, green, purple        %fp: red, orange, /                %fn: gray, black, /
    colors = [0.30 0.75 0.93; 0 0.5 0; 0.5 0 0.5; 0.9 0 0; 0.929 0.694 0.125; 0 0 0; 0.5 0.5 0.5; 0 0 0; 0 0 0];
    colors = mat2cell(colors,ones(size(colors,1),1),3);
    stackData = []; 

    current_sub = data.(fn{is});
    stackData(:,1,:) = current_sub.all_ref';
    stackData(:,2,:) = current_sub.all_ref_in'; 
    stackData(:,3,:) = current_sub.all_KS'; 
    h = F4_plotBarStackGroups(stackData, groupLabels,subject, is); hold on
    width = 0.3;
    set(h,{'FaceColor'},colors)
    legend('ref-correct','ref-incorrect','ref-unmatched','ref-correct at \Delta z = 10\mum','ref-incorrect at \Delta z = 10\mum','ref-unmatched at \Delta z = 10\mum','non-ref at \Delta z = 10\mum')
    for ir = 1:day
        switch subject{is}
            case 'AL031'
                groupBins = [1 1/3+3 6 7];
                xpos = groupBins(ir)-width-width/2;
                xpos_in = groupBins(ir)-width/2;
            case 'AL032'
                groupBins = [1 2/3+3 4.66 8];
                xpos = groupBins(ir)-width-width/2;
                xpos_in = groupBins(ir)-width/2;
            case 'AL036'
                xpos = ir-width-width/2;
                xpos_in = ir-width/2;
        end
        text(xpos,current_sub.corr_ref(1,ir)-4,num2str(round(current_sub.rate_ref(1,ir),2)),"FontSize",16); hold on
        text(xpos_in,current_sub.corr_ref_in(1,ir)-4,num2str(round(current_sub.rate_ref_in(1,ir),2)),"FontSize",16); hold on
    end
    ax = gca;
    ax.FontSize = 16;
    set(ax,'XLim',[0 9])
    set(ax,'TickDir','out');
    xlabel('Days matched','FontSize',18,'FontWeight','Bold','Color','k')
    ylabel('Number of units','FontSize',18,'FontWeight','Bold','Color','k')
    title(sprintf('%s Accuracy',subject{is}))
end






