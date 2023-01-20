clear all
addpath('D:\Data\Pipeline\npy\npy-matlab\npy-matlab'); %npy path
addpath('C:\Users\labadmin\Desktop\Neuron_Tracking-main\wave_metrics');
addpath(genpath('C:\Users\labadmin\Desktop\Neuron Tracking Pipeline'));

% ----------------- User input: path -----------------
input_struct.rootpath = 'Z:\'; %all data
input_struct.subject = 'AL031';
input_struct.sorting_path = ['Z:\',input_struct.subject,'_out\results'];
input_struct.EMD_input_dir = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\EMD_input'; %stores EMD input pre-correction
input_struct.ref_path = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\reference'; %path to save reference results
input_struct.clu_data = fullfile(input_struct.ref_path,[input_struct.subject,'_vfp.mat']); %reference data path
input_struct.psth_data = fullfile(input_struct.ref_path,[input_struct.subject,'_PSTH.mat']);
input_struct.output_dir = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\output';

% ----------------- User input: setting -----------------
input_struct.day = 5; %number of days to be compared
input_struct.shank = 1; %number of shanks
if input_struct.day < 2
    error('Need at least 2 days.');
elseif isinteger(input_struct.day)
    error('Please enter integer days.');
end
input_struct.plot_mask = logical([0,0,0,0,0]);
%input_struct.plot_mask = logical([1,1,1,1,1]); %plot order: plot_ref_hist, plot_ref_match, plot_KSgood_hist, plot_KSgood_ref_hist, plot_KSgood_match_hist

% ----------------- Set to default, can change if necessary -----------------
input_struct.n = 8; %ref parameter: sim_score range
input_struct.max_dist = 36; %ref parameter: um, 8 possible sites
%input_struct.dim_mask = logical([1,1,1,0,0,0,0,0,0,0]);   %for just position
input_struct.dim_mask = logical([1,1,1,0,0,0,0,0,0,1]);   %for positions + waveform l2
%input_struct.dim_mask = logical([0,0,0,0,0,0,0,0,0,1]);    %for just wf
input_struct.dim_mask_corrected = logical([1,1,1,0,0,0,0,0,0,1]);   %for positions + waveform l2
%input_struct.dim_mask_corrected = logical([1,1,1,0,0,0,0,0,0,0]);
%input_struct.dim_mask_corrected = logical([0,0,0,0,0,0,0,0,0,1]);
input_struct.location = 'wf'; %'wf'; %'Boussard';
input_struct.mode_est_method = 'single'; %'single'; 'merge'
input_struct.l2_weights = 0;
all_l2_weights = [1500];

%input_struct.corr_method = 'direct'; %correction by 'row' or 'direct'
%input_struct.ref_criterion = 'KS_good'; %'amp'
%input_struct.drift_est_method = 'kernel'; %'kernel'; %'mode';, 'mode_allShank', need to remove





%% Pre correction
% % calculate cluster centroids if haven't
% calcUnitCentroids(rootpath,day,shank);
for iw = 1:length(all_l2_weights)

    input_struct.l2_weights = all_l2_weights(iw);

    pair_output.z_mode = [];
    output_struct.KS = [];
    diff_comb = [];
    % % Plot reference pairs change: KS, KW test, ref
    % color = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
    % figure()
    % xlabel('Day')
    % ylabel('Number of clusters')
    % plot(output_struct.KS(:,1),'-.','LineWidth',1,'Color',color(3,:)); hold on;
    % plot(output_struct.KW(:,1),'--','LineWidth',1,'Color',color(2,:)); hold on;
    % plot([2 3 4],output_struct.ref_pre(:,1),'-p','LineWidth',2,'Color',color(1,:)); hold on;
    % set(gca,'xtick',[1:1:4],'XTickLabel',{'day1','day2','day3','day4'})
    % legend('KS','KW sig','ref pairs')
    % legend('Location','best')
    % for ir = 1:size(output_struct.KW,1)
    %     text(ir,output_struct.KS(ir,1)+1,num2str(output_struct.KS(ir,1)),"FontSize",13);
    %     text(ir,output_struct.KW(ir,1)+1,num2str(output_struct.KW(ir,1)),"FontSize",13);
    % end
    % for ir = 1:size(output_struct.ref_pre,1)
    %     text(ir+1,output_struct.ref_pre(ir,1)+3,num2str(output_struct.ref_pre(ir,1)),"FontSize",13);
    % end
    % title('Example change of clusters, shank0')


    %
    if ~exist(fullfile(input_struct.ref_path, [input_struct.subject,'_stimulus_times.mat']),'file')
        stim_time_data(input_struct);
        calculate_PSTH(input_struct);
        vfp_KW(input_struct);
    end
    input_struct.day_label = load(fullfile(input_struct.ref_path, [input_struct.subject,'_stimulus_times.mat'])).day_label;

    % loop all days and shanks
    d = dir(input_struct.sorting_path); %list all files in directory
    d = d([d.isdir]); %list sub-directories only
    d = d(~ismember({d.name},{'.', '..'})); %exclude . and .. directories
    for id = 1:input_struct.day-1
        date = input_struct.day_label(id+1);
        idx = find(contains({d.name},date));
        fday1 = fullfile(input_struct.sorting_path, d(find(contains({d.name},input_struct.day_label(1)))).name); % all compare to d1
        fday2 = fullfile(input_struct.sorting_path, d(idx).name);
        dday1 = dir(fday1);
        dday1 = dday1([dday1.isdir]);
        dday1 = dday1(~ismember({dday1.name},{'.', '..'}));
        dday2 = dir(fday2);
        dday2 = dday2([dday2.isdir]);
        dday2 = dday2(~ismember({dday2.name},{'.', '..'}));
        for ish = 1:input_struct.shank %4 shanks
            % Find reference set
            [output_struct,pair_output] = find_reference(input_struct,'pre',output_struct,pair_output,id,ish);
            % plot reference pairs
            plot_ref_hist(input_struct, pair_output,'pre', id, ish);

            fname1 = fullfile(fday1, dday1(ish).name);
            fname2 = fullfile(fday2, dday2(ish).name);
            phydir1 = [fname1,'\imec',num2str(ish-1),'_ks25'];
            locdir1 = [fname1,'\Localization\position_results_files_merged'];
            phydir2 = [fname2,'\imec',num2str(ish-1),'_ks25']; %the later day
            locdir2 = [fname2,'\Localization\position_results_files_merged'];
            %EMD_input_name = ['d',num2str(id+1),'d',num2str(id),'_sh',num2str(ish-1),'_locwf_emd_input.mat'];
            EMD_input_name = ['d',num2str(id+1),'d1_sh',num2str(ish-1),'_locwf_emd_input.mat']; %d1_based

            % create input file for EMD correction
            EMD_input(input_struct, output_struct, pair_output, phydir1, phydir2, locdir1, locdir2, EMD_input_name,'pre',ish, id);
            % plot reference pair and KSgood units
            plot_ref_match(input_struct,EMD_input_name,'pre',id,ish);

            % Run EMD to find correction amount
            %matched_filename = ['d',num2str(id+1),'d',num2str(id),'_sh',num2str(ish-1),'_locwf_emd_output.mat'];
            matched_filename = ['d',num2str(id+1),'d1_sh',num2str(ish-1),'_locwf_emd_output.mat'];
            pair_output = EMD_match(input_struct,pair_output,EMD_input_name,matched_filename,'pre',id,ish);
            plot_ref_match(input_struct, matched_filename,'pre',id,ish); %plot matched

            % estimate and plot histograms
            [diffZ,edges] = z_estimate(input_struct.EMD_input_dir, matched_filename);
            plot_KSgood_hist(input_struct,matched_filename,'pre',id,ish);
            plot_KSgood_ref_hist(input_struct,matched_filename,'pre',id,ish);
            plot_KSgood_match_hist(input_struct,matched_filename,'pre',id,ish);


            %single-day
            switch input_struct.mode_est_method
                case 'single'
                    pair_output.z_mode = kernelModeEstimate(diffZ);

                    % Correct drift for centroid locations
                    switch input_struct.location
                        case 'Boussard'
                            correct_drift(pair_output.z_mode,locdir2,id,ish);
                    end
            end


            % Post-correction
            % Correct drift for reference data
            [output_struct, pair_output] = find_reference(input_struct,'post',output_struct,pair_output,id,ish);
            % plot reference pairs
            plot_ref_hist(input_struct,pair_output, 'post', id, ish);

            EMD_input_name_corrected = ['d',num2str(id+1),'d1_sh',num2str(ish-1),'_locwf_emd_input_corrected.mat']; %d1_based

            % create input file for EMD correction
            EMD_input(input_struct, output_struct, pair_output, phydir1, phydir2, locdir1, locdir2, EMD_input_name_corrected,'post',ish, id);
            % plot reference pair and KSgood units
            plot_ref_match(input_struct,EMD_input_name_corrected,'post',id,ish);

            matched_filename_corrected = ['d',num2str(id+1),'d1_sh',num2str(ish-1),'_locwf_emd_output_corrected.mat'];
            pair_output = EMD_match(input_struct,pair_output,EMD_input_name_corrected,matched_filename_corrected,'post',id,ish);
            plot_ref_match(input_struct, matched_filename,'pre',id,ish); %plot matched

            [diffZ_post,edges_post] = z_estimate(input_struct.EMD_input_dir, matched_filename_corrected);

            % Find z-correction amount
            pair_output.diffZ_post = diffZ_post;
            pair_output.z_mode_post = kernelModeEstimate(diffZ_post);

            % save pair_output
            if ~exist(input_struct.output_dir, 'dir')
                mkdir(input_struct.output_dir);
            end

            pairout_name = [input_struct.subject,'d1',num2str(id+1),'sh',num2str(ish),'_pair_output_',num2str(all_l2_weights(iw)),'.mat'];
            save(fullfile(input_struct.output_dir, pairout_name),'pair_output');


            %         % figure 3: drift plot
            %         plotFlow(input_struct.EMD_input_dir, matched_filename,pair_output.z_mode,'pre',id, ish);


            %         % combine histograms
            %         if ish ~= 2 % NEED TO CHANGE
            %             diff_comb = [diff_comb; diffZ]; %combine all diffZ
            %         end
            %     end
            %
            %     %merged
            %     switch input_struct.mode_est_method
            %         case 'merge'
            %             diffZ = diff_comb;
            %             pair_output.z_mode = kernelModeEstimate(diffZ);
            %
            %             for ish = 1:input_struct.shank %4 shanks
            %                 fname1 = fullfile(fday1, dday1(ish).name);
            %                 fname2 = fullfile(fday2, dday2(ish).name);
            %                 phydir1 = [fname1,'\imec',num2str(ish-1),'_ks25'];
            %                 locdir1 = [fname1,'\Localization\position_results_files_merged'];
            %                 phydir2 = [fname2,'\imec',num2str(ish-1),'_ks25']; %the later day
            %                 locdir2 = [fname2,'\Localization\position_results_files_merged'];
            %
            %                 % Correct drift for centroid locations
            %                         correct_drift(pair_output.z_mode,locdir2, id, ish);
            %             end
        end
    end

    input_name = [input_struct.subject,'input_setting_',num2str(all_l2_weights(iw)),'.mat'];
    save(fullfile(input_struct.output_dir, input_name),"input_struct");
    output_name = [input_struct.subject,'all_output_',num2str(all_l2_weights(iw)),'.mat'];
    save(fullfile(input_struct.output_dir, output_name),"output_struct");

end
% plot accuracy line
% figure()
% pair_good = output_struct.corr_pre + output_struct.fp_pre + output_struct.fn_pre;
% color = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
% xlabel('Comparison days','FontSize',18)
% ylabel('Accuracy','FontSize',18)
% set(gca,'xtick',[1:1:3],'XTickLabel',{'day2-1','day3-1','day4-1'})
% title('Pre-correction correct pairs')
% for ip = 1:size(output_struct.corr_pre,2)
%     plot(output_struct.corr_pre(:,ip),'-p','LineWidth',2,'Color',color(ip,:)); hold on;
%     plot(pair_good(:,ip),'--','LineWidth',1,'Color',color(ip,:)); hold on;
% end
% legend('shank 0','shank 0 all','shank 1','shank 1 all','shank 2','shank 2 all','shank 3','shank 3 all')
% for ir = 1:size(output_struct.pair_results,1)
%     for ic = 1:size(output_struct.pair_results,2)
%         numPair(ir,ic) = length(output_struct.pair_results{ir,ic});
%         if ic == 2 | ic == 4
%             text(ir,output_struct.corr_pre(ir,ic),[num2str(numPair(ir,ic)*output_struct.corr_pre(ir,ic)),'/',num2str(numPair(ir,ic))],"FontSize",15);
%         else
%             text(ir,output_struct.corr_pre(ir,ic)-0.02,[num2str(numPair(ir,ic)*output_struct.corr_pre(ir,ic)),'/',num2str(numPair(ir,ic))],"FontSize",15);
%         end
%     end
% end




%% plot

% % plot bar results
% color = [0.4660 0.6740 0.1880; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.8 0.8 0.8];
% all_Ref = pair_output.ref_pre;
% for ip = 1:input_struct.shank
%     figure()
%     result = [output_struct.corr_pre(:,ip)'; output_struct.fp_pre(:,ip)'; output_struct.fn_pre(:,ip)']';
%     all_good = sum(result,2);
%     result(:,4) = all_Ref(:,ip)-all_good;
%     ba = bar(result,'stacked','FaceColor','flat'); hold on
%     ba(1).CData = color(1,:);
%     ba(2).CData = color(2,:);
%     ba(3).CData = color(3,:);
%     ba(4).CData = color(4,:);
%     %result_added = [result(:,1)'; (result(:,1)+result(:,2))'; (result(:,1)+result(:,2)+result(:,3))']';
%     plot(all_Ref(:,ip),'LineWidth',2,'Color',[0.1 0.1 0.1]); hold on
%     plot(all_good,'LineWidth',2,'Color',[0.1 0.1 0.1]); hold on
%     plot(output_struct.corr_pre(:,ip),'LineWidth',2,'Color',[0.1 0.1 0.1])
%     title(sprintf('Pre-correction results sh%d', ip))
%     set(gca,'xtick',[1:1:3],'XTickLabel',{'day2-1','day3-1','day4-1'})
%     legend('correct','FP','FN','MUA')
%     for ir = 1:size(output_struct.pair_results,1)
%         text(ir,all_Ref(ir,ip)+1,num2str(all_Ref(ir,ip)),"FontSize",13);
%         text(ir,all_good(ir)+1,num2str(all_good(ir)),"FontSize",13);
%         text(ir,output_struct.corr_pre(ir,ip)+1,num2str(output_struct.corr_pre(ir,ip)),"FontSize",13);
%         text(ir,1,num2str(round(output_struct.corr_pre(ir,ip)/all_good(ir),2)),"FontSize",13);
%     end
%     hold off
% end
%
% % plot bar results
% color_post = [0 1 0; 1 0 0; 1 1 0; 0.8 0.8 0.8];
% all_Ref_post = pair_output.ref_post;
% for ip = 1:input_struct.shank
%     figure()
%     result_post = [output_struct.corr_post(:,ip)'; output_struct.fp_post(:,ip)'; output_struct.fn_post(:,ip)']';
%     all_good = sum(result_post,2);
%     result_post(:,4) = all_Ref_post(:,ip)-all_good;
%     ba = bar(result_post,'stacked','FaceColor','flat'); hold on
%     ba(1).CData = color_post(1,:);
%     ba(2).CData = color_post(2,:);
%     ba(3).CData = color_post(3,:);
%     ba(4).CData = color_post(4,:);
%     plot(all_Ref_post(:,ip),'LineWidth',1,'Color',[0.1 0.1 0.1]); hold on
%     plot(all_good,'LineWidth',1,'Color',[0.1 0.1 0.1])
%     plot(output_struct.corr_post(:,ip),'LineWidth',1,'Color',[0.1 0.1 0.1])
%     title(sprintf('Post-correction results sh%d', ip))
%     set(gca,'xtick',[1:1:3],'XTickLabel',{'day2-1','day3-1','day4-1'})
%     legend('correct','FP','FN','MUA')
%     for ir = 1:size(output_struct.pair_results_post,1)
%         text(ir,all_Ref_post(ir,ip)+1,num2str(all_Ref_post(ir,ip)),"FontSize",13);
%         text(ir,all_good(ir)+1,num2str(all_good(ir)),"FontSize",13);
%         text(ir,output_struct.corr_post(ir,ip)+1,num2str(output_struct.corr_post(ir,ip)),"FontSize",13);
%         text(ir,1,num2str(round(output_struct.corr_post(ir,ip)/all_good(ir),2)),"FontSize",13);
%     end
%     hold off
% end







% plot(output_struct.corr_post,'-p','LineWidth',2)
% legend('shank 0','shank 1','shank 2','shank 3')
% xlabel('Comparison days','FontSize',18)
% ylabel('Accuracy','FontSize',18)
% ylim([0 1])
% yticks([0:0.1:1])

% Bar plot of results
% figure()
% result = [output_struct.corr_post(:,1)' output_struct.fp_post(:,1)' output_struct.fn_post(:,1)'];
% bar(result,'stacked')
% title('Post-correction results')
% set(gca,'xtick',[1:1:3],'XTickLabel',{'day2-1','day3-1','day4-1'})
% for ir = 1:size(output_struct.pair_results_post,1)
%     for ic = 1:size(output_struct.pair_results_post,2)
%         nPair = length(output_struct.pair_results_post{ir,ic});
%         if ic == 2 | ic == 4
%             text(ir,output_struct.corr_post(ir,ic),[num2str(nPair*output_struct.corr_post(ir,ic)),'/',num2str(nPair)],"FontSize",15);
%         else
%             text(ir,output_struct.corr_post(ir,ic)-0.02,[num2str(nPair*output_struct.corr_post(ir,ic)),'/',num2str(nPair)],"FontSize",15);
%         end
%     end
% end

