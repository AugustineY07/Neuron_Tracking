function plot_ref_match(input_struct,input_name,stage,id,ish)
% function to plot unit locations for two days in XZ, colored by one other
% metric (hard coded to amplitude, see line 10.
% The 'flow' lines are the reference matches for an input file for
% EMD_main, and the emd-found pairs for an output file from EMD_main.
% input argument is th full path to the .mat file

% if isempty(varargin)
%     [dat_name, dat_folder] = uigetfile('*.mat','Select file of pairs and matches');
%     dat_fullpath = fullfile(dat_folder,dat_name);
% else
%     inputCell = varargin(1);
%     dat_fullpath = inputCell{1};
% end

opt = input_struct.plot_mask(2);
EMD_input_dir = input_struct.EMD_input_dir;
dat = load(fullfile(EMD_input_dir,input_name));

% scatter plot colored by amplitude
switch opt
    case 1
        h1 = figure();
        h1.Name = sprintf('%s-correction flow day %d-1 shank %d',stage,id+1,ish);
        h1.Position = [110,40,310,1110];

        color_ind = 4;
        c = dat.f1(:,color_ind);
        scatter(dat.f1(:,1), dat.f1(:,2),50,c,'o');
        hold on;
        c = dat.f2(:,color_ind);
        scatter(dat.f2(:,1), dat.f2(:,2),50,c,'square','filled');

        if isfield(dat,'f2_same_ind')
            match = dat.f2_same_ind;
        else
            match = dat.f2_emd_ind;
        end

        % loop over pairs, drawing lines
        nMatch = sum(~isnan(match))
        for i = 1:length(dat.f1)
            cm = match(i);
            if ~isnan(cm)
                plot([dat.f2(cm,1),dat.f1(i,1)],[dat.f2(cm,2),dat.f1(i,2)],"black");
            end
        end
        xlim([-20,50]);
        ylim([2800,3600]);
        colorbar
        caxis([80,300])
        hold off;
end







