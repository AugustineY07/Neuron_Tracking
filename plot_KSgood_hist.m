function plot_KSgood_hist(input_struct,matched_filename,stage,id,ish)

opt = input_struct.plot_mask(3);
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

        % all KSgood histograms
        h3 = figure();
        h3.Name = sprintf('%s-correction day %d-1 shank %d zdiff_his',stage,id+1,ish);
        h3.Units = 'centimeters';
        h3.Position = [8,6,10,8];

        histogram(diffZ,edges);
        ax = gca;
        ax.FontSize = 18;
        xlim([-80,80])
        xlabel('EMD matched Vertical-shift (\mum)','FontSize',18,'FontWeight','bold','Color','k')
        ylabel('Number of matches','FontSize',18,'FontWeight','bold','Color','k')
end