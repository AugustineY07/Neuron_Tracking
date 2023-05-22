function vfp_KW(input_struct)
% Apply Kruskal-Wallis test to select visual fingerprints of a single day that are significantly different from each other

day = input_struct.day;
shank = input_struct.shank;
input_path = input_struct.rf_path;
output_path = input_struct.rf_path;
input_name = input_struct.psth_data;
start_day = input_struct.start_day;
compare = input_struct.compare;
output_name = input_struct.clu_data;

load(fullfile(input_path,input_name));

% switch compare
%     case '1-based'
%         output_name = [input_struct.subject, '_vfp.mat'];
%     otherwise
%         output_name = [input_struct.subject, '_vfp_skip3_', compare, '.mat'];
% end

% compute KW p-value for all fingerprints
for ish = 1:size(rsp,2)
    for iday = 1:size(rsp,1)
        n_good = length(good_clu(iday,ish).clu);
        for i = 1:n_good
            sample = [];
            tcount = 0;
            for itrial = 1:size(rspByTrial,3)
                if isempty(rspByTrial{iday,ish,itrial})
                    continue
                else
                    tcount = tcount + 1;
                    sample(tcount,:) = rspByTrial{iday,ish,itrial}(i,:);
                    p{iday,ish}(i) = kruskalwallis(sample,[],'off'); %p value for every cluster in every day
                end
            end
        end
        p_sig{iday,ish} = find(p{iday,ish} < 0.01); %signifiant col idx among input clusters 
        clu_sig{iday,ish} = good_clu(iday,ish).clu(p_sig{iday,ish}); %cluster idx of significant clusters
        % rsp and PSTH of KW significant clusters
        rsp_sig{iday,ish} = rsp{iday,ish}(p_sig{iday,ish},:);
        PSTH_sig{iday,ish} = PSTHsmoothed{iday,ish}(p_sig{iday,ish},:);
    end
end



% clusters passes KW test
for ish = 1:shank %loop on shanks  
    for id = 1:day-1 % loop on days

        shank = [];
        shank2 = [];
        shank(:,1) = good_clu(id,ish).clu(p_sig{id,ish}); %current day, significant cluster labels
        shank(:,2) = good_clu(id,ish).ch(p_sig{id,ish});
        shank(:,3) = good_clu(id,ish).x(p_sig{id,ish});
        shank(:,4) = good_clu(id,ish).y(p_sig{id,ish});
        shank(:,5) = good_clu(id,ish).amp(p_sig{id,ish});
        shank(:,6) = [1:1:length(good_clu(id,ish).clu(p_sig{id,ish}))]'; %original order in rsp and PSTH (if rank by clu label)
        shank(:,7) = good_clu(id,ish).ISI(p_sig{id,ish});
        shank(:,8) = p{id,ish}(p_sig{id,ish})';
        shank = sortrows(shank,[4 3]); %sort by y location, then by x location

        shank2(:,1) = good_clu(id+1,ish).clu(p_sig{id+1,ish}); %next day
        shank2(:,2) = good_clu(id+1,ish).ch(p_sig{id+1,ish});
        shank2(:,3) = good_clu(id+1,ish).x(p_sig{id+1,ish});
        shank2(:,4) = good_clu(id+1,ish).y(p_sig{id+1,ish});
        shank2(:,5) = good_clu(id+1,ish).amp(p_sig{id+1,ish});
        shank2(:,6) = [1:1:length(good_clu(id+1,ish).clu(p_sig{id+1,ish}))]';
        shank2(:,7) = good_clu(id+1,ish).ISI(p_sig{id+1,ish});
        shank2(:,8) = p{id+1,ish}(p_sig{id+1,ish})';
        shank2 = sortrows(shank2,[4 3]);

        %for save
        good_clu_KW{id,ish} = shank;
        good_clu_KW{id+1,ish} = shank2;
    end
end



% calculate rsp and PSTH, need to choose START day
count = 0;
for ish = 1:size(rsp,2) %loop on shanks
    for id = start_day+1:day %loop on other day
        % calculate correlation of rsp and PSTH between start day and later days
        rsp_sort = rsp_sig{start_day,ish}(good_clu_KW{start_day,ish}(:,6),:); %rearrange rsp at start day based on physical location. now the same order as in good_clu_KW
        PSTH_sort = PSTH_sig{start_day,ish}(good_clu_KW{start_day,ish}(:,6),:);
        rsp_sort2 = rsp_sig{id,ish}(good_clu_KW{id,ish}(:,6),:); %rearrange later days based on physical location
        PSTH_sort2 = PSTH_sig{id,ish}(good_clu_KW{id,ish}(:,6),:);
        sim_rsp_sig{id,ish} = corr(rsp_sort', rsp_sort2');
        sim_PSTH_sig{id,ish} = corr(PSTH_sort', PSTH_sort2');
        simScore_sig{id,ish} = sim_rsp_sig{id,ish} + sim_PSTH_sig{id,ish}; %contains negative corr and NaN

        % calculate correlation for all clusters regardless of KW significance
        sim_rsp{id,ish} = corr(rsp{start_day,ish}', rsp{id,ish}');
        sim_PSTH{id,ish} = corr(PSTHsmoothed{start_day,ish}', PSTHsmoothed{id,ish}');
        simScore{id,ish} = sim_rsp{id,ish} + sim_PSTH{id,ish};

        switch compare
            case '1-based'
                % calculate correlation of rsp and PSTH between day1 and other days
                rsp_sort = rsp_sig{1,ish}(good_clu_KW{1,ish}(:,6),:); %rearrange d1 based on physical location
                PSTH_sort = PSTH_sig{1,ish}(good_clu_KW{1,ish}(:,6),:);
                rsp_sort2 = rsp_sig{id+1,ish}(good_clu_KW{id+1,ish}(:,6),:); %rearrange others based on physical location
                PSTH_sort2 = PSTH_sig{id+1,ish}(good_clu_KW{id+1,ish}(:,6),:);
                sim_rsp_sig{id,ish} = corr(rsp_sort', rsp_sort2');
                sim_PSTH_sig{id,ish} = corr(PSTH_sort', PSTH_sort2');
                simScore_sig{id,ish} = sim_rsp_sig{id,ish} + sim_PSTH_sig{id,ish}; %contains negative corr and NaN
            case 'next'
                % calculate correlation of rsp spike count and mean PSTH between the same shank different day
                rsp_sort = rsp_sig{id,ish}(good_clu_KW{id,ish}(:,6),:); %rearrange based on physical location
                PSTH_sort = PSTH_sig{id,ish}(good_clu_KW{id,ish}(:,6),:);
                rsp_sort2 = rsp_sig{id+1,ish}(good_clu_KW{id+1,ish}(:,6),:); %rearrange based on physical location
                PSTH_sort2 = PSTH_sig{id+1,ish}(good_clu_KW{id+1,ish}(:,6),:);
                sim_rsp_sig{id,ish} = corr(rsp_sort', rsp_sort2');
                sim_PSTH_sig{id,ish} = corr(PSTH_sort', PSTH_sort2');
                simScore_sig{id,ish} = sim_rsp_sig{id,ish} + sim_PSTH_sig{id,ish}; %contains negative corr and NaN


        end


%         % correlation heatmap
%         count = count+1;
%         figure(count)
%         h = heatmap(simScore_sig{id,ish});
%         h.FontSize = 12;
%         h.XLabel = ['Day ',sprintf('%d',id+1),' Shank ', sprintf('%d', ish)];
%         h.YLabel = ['Day ',sprintf('%d',id),' Shank ', sprintf('%d', ish)];
%         h.Title = 'SimScore';
%         h.Colormap = parula;
% 
%         xPos = [strcat(string(shank2(:,3)), string(shank2(:,4)))];
%         yPos = [strcat(string(shank(:,3)), string(shank(:,4)))];
%         h.XDisplayLabels = xPos;
%         h.YDisplayLabels = yPos;
    end
end

p_value = p; 
save(fullfile(output_path,output_name), 'good_clu_KW', 'p_sig','clu_sig', 'rsp_sig', 'sim_rsp_sig', 'PSTH_sig', 'sim_PSTH_sig','simScore_sig','sim_rsp','sim_PSTH','simScore','p_value')
