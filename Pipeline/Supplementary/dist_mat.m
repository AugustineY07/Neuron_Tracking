%Plot the EMD cost, location cost and waveform cost distance matrix
%Perform (optional) Mann-Whitney U Test between datasets
function dist_mat(all_output,numData)
cost_wf_all = zeros(numData,numData);
cost_loc_all = zeros(numData,numData);
cost_all = zeros(numData,numData);

shank = 1;

for id1 = 1:numData-1 %first day
    for id2 = id1+1:numData %second day
        cost_wf = 0;
        cost_loc = 0;
        cost = 0;
        pair_output = all_output(id1).output;
        C_wf = pair_output.C_wf_post(:); %vectorize distance matrix
        C_loc = pair_output.C_physical_post(:);
        num1 = pair_output.KSgood_f1;
        num2 = pair_output.KSgood_f2;
        KS_norm = min([num1,num2]);

        % cost
        cost_wf = cost_wf + sum(C_wf.*pair_output.x_post) / (KS_norm*shank); %normalized cost by KSgood unit number
        cost_loc = cost_loc + sum(C_loc.*pair_output.x_post) / (KS_norm*shank);
        cost = cost + pair_output.cost_post / (KS_norm*shank);

    end
    cost_wf_all(id1,id2) = cost_wf;
    cost_loc_all(id1,id2) = cost_loc;
    cost_all(id1,id2) = cost;
end

% plot heatmap
h = figure();
ax = gca;
ax.FontSize = 18; %tick font
ax.FontWeight = 'Bold';

subplot(1,3,1)
heatmap(round(cost_wf_all))
xlabel('Second dataset index')
ylabel('First dataset index')
title('Waveform cost in each comparison')

subplot(1,3,2)
heatmap(round(cost_loc_all))
xlabel('Second dataset index')
ylabel('First dataset index')
title('Location cost in each comparison')

subplot(1,3,3)
heatmap(round(cost_all))
xlabel('Second dataset index')
ylabel('First dataset index')
title('Cost in each comparison')


end