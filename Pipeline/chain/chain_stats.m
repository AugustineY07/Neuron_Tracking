% calculate L2, FR, and locations of all chains
function [L2_weight,fr_all,fr_change,x_loc,z_loc] = chain_stats(all_input,all_output,full_chain,numData)

for ichain = 1:size(full_chain,1)
    for id = 1:numData-1
        % current clu
        clu_label1 = full_chain(ichain,id);
        clu_label2 = full_chain(ichain,id+1);

        % Calcuate L2 weight between clu_label1,clu_label2
        L2 = reshape(all_output(id).output.L2, all_output(id).output.KSgood_f1, all_output(id).output.KSgood_f2);
        L2_wt = L2'.*all_output(id).output.P_post;

        pair_results = load(fullfile(all_input(id).input.EMD_path,['EMD_post',num2str(id),'.mat']));
        clusters1 = pair_results.f1_labels; %get cluster label
        clusters2 = pair_results.f2_labels;
        idx1 = find(clusters1 == clu_label1); %get cluster idx in all KSgood
        idx2 = find(clusters2 == clu_label2);
        L2_weight(ichain,id) = L2_wt(idx2,idx1); %get L2 weight





        % Calculate firing rate change
        % read excel summary
        tb1 = readtable(fullfile(all_input(id).input.input_path, all_input(id).input.data_path1, 'metrics.csv'));
        tb2 = readtable(fullfile(all_input(id).input.input_path, all_input(id).input.data_path2, 'metrics.csv'));

        % find FR
        fr{id} = table2array(tb1(:,2)); %get unit firing rate
        fr_ave(id) = sum(fr{id})/size(tb1,1); %average across num clu
        fr{id+1} = table2array(tb2(:,2));
        fr_ave(id+1) = sum(fr{id+1})/size(tb2,1);


        % find all clusters
        clu1 = table2array(tb1(:,1))+1;
        clu2 = table2array(tb2(:,1))+1;

        % calculate fold change
        idx_row_ref = find(clu_label1 == clu1); %idx in all day 1 clusters
        fr_all(ichain,id) = fr{id}(idx_row_ref);
        fr_fold_all(ichain,id) = fr{id}(idx_row_ref)/fr_ave(id); %how is unit is firing compare to average of this day
        idx_row_ref2 = find(clu_label2 == clu2);
        fr_all(ichain,id+1) = fr{id+1}(idx_row_ref2);
        fr_fold_all(ichain,id+1) = fr{id+1}(idx_row_ref2)/fr_ave(id+1);
        fr_change(ichain,id) = (fr_fold_all(ichain,id+1)-fr_fold_all(ichain,id))/fr_fold_all(ichain,id); %percentage change of FR compare to prior day




        % find z locations
        z_drift(ichain,id+1) = all_output(id).output.z_mode;
        pair_results = load(fullfile(all_input(id).input.EMD_path, all_input(id).input.filename_post));
        f1 = pair_results.f1; %get all location data
        f2 = pair_results.f2;
        f1_labels = pair_results.f1_labels;
        f2_labels = pair_results.f2_labels;
        idx1 = find(f1_labels == clu_label1); %find clu index in all units
        idx2 = find(f2_labels == clu_label2);
        z_loc(ichain,id) = f1(idx1,2); %all z locations, col = days
        x_loc(ichain,id) = f1(idx1,1);
        if id == numData-1
            z_loc(ichain,id+1) = f2(idx2,2); 
            x_loc(ichain,id+1) = f2(idx2,1);
        end
    end
end

end