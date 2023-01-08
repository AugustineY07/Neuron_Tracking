function plot_KSgood_ref_hist(input_struct,matched_filename,stage,id,ish)

opt = input_struct.plot_mask(4);
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

        % KSgood reference pair histogram
        f1_KSgood = find(ismember(f1_labels,pair_results(:,1))); %KSgood row in f1
        f2_KSgood(:,1) = find(ismember(f2_labels,pair_results(:,2))); %KSgood row in f2
        f2_KSgood(:,2) = f2_labels(f2_KSgood(:,1)); %KSgood label
        [~,f1order] = ismember(pair_results(:,2),f2_KSgood(:,2));
        f2_KSgood_sort = f2_KSgood(f1order,:);
        diffZ_good = f2(f2_KSgood_sort(:,1),2) - f1(f1_KSgood,2);

        h4 = figure();
        histogram(diffZ_good,edges);
        ax = gca;
        ax.FontSize = 18;
        xlim([-80,80])
        title(sprintf('%s-correction day %d-1 shank %d KSgood+ref',stage,id+1,ish))
        xlabel('EMD matched Vertical-shift (\mum)','FontSize',18,'FontWeight','bold','Color','k')
        ylabel('Number of matches','FontSize',18,'FontWeight','bold','Color','k')
end
