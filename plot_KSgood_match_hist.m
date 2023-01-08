function plot_KSgood_match_hist(input_struct,matched_filename,stage,id,ish)

opt = input_struct.plot_mask(5);
inputDir = input_struct.EMD_input_dir;
fileName = matched_filename;

switch opt
    case 1
        load(fullfile(inputDir,fileName));
        hasMatch = ~isnan(f2_emd_ind);
        f1_matched = f1(hasMatch,:);
        f2_matched = f2(f2_emd_ind(hasMatch),:);
        diffZ = f2_matched(:,2) - f1_matched(:,2);
        binsize = 4;
        edges = (-100:binsize:100);

        % KSgood+matching reults histogram
        idx_used = [];
        for iks = 1:size(pair_results,1)
            if pair_results(iks,4) == -1
                f2_KSmatch(iks,:) = 0;
            elseif pair_results(iks,5) == -1
                f2_KSmatch(iks,:) = 0;
            else
                f1_KSmatch(iks,1) = find(ismember(f1_labels,pair_results(iks,5))); %matched row in f1
                f2_KSmatch(iks,1) = find(ismember(f2_labels,pair_results(iks,4))); %KSmatched row in f2
                f2_KSmatch(iks,2) = f2_labels(f2_KSmatch(iks,1)); %KSgood label
                idx_used = [idx_used iks];
            end
        end
        diffZ_KSmatch = f2(f2_KSmatch(idx_used,1),2) - f1(f1_KSmatch(idx_used),2);

        figure()
        histogram(diffZ_KSmatch,edges);
        ax = gca;
        ax.FontSize = 18;
        xlim([-80,80])
        title(sprintf('%s-correction day %d-1 shank %d KSgood+matched',stage,id+1,ish))
        xlabel('EMD matched Vertical-shift (\mum)','FontSize',18,'FontWeight','bold','Color','k')
        ylabel('Number of matches','FontSize',18,'FontWeight','bold','Color','k')
end

