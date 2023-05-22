function stim_time_data(input_struct)
% This file reads and organizes stimulus id and times
% output: sorted_response.mat

input_path = fullfile(input_struct.rootpath,input_struct.subject);
output_path = input_struct.ref_path;
out_name = [input_struct.subject,'_stimulus_times.mat'];
day = input_struct.day;

% read in stimulus data
d = dir(input_path); %list all files in directory
d = d([d.isdir]); %list sub-directories only
d = d(~ismember({d.name},{'.', '..'})); %exclude . and .. directories
icount = 1;
iloop = 1;
while icount <= day && iloop <= length(d) %number of days to compare
    d1 = fullfile(input_path, d(iloop).name,'alignments');
    sid = fullfile(d1, ['stimulus_ids_', input_struct.subject, '_', d(iloop).name, '.mat']);
    stim = fullfile(d1, ['stimulus_times_', input_struct.subject, '_', d(iloop).name, '.mat']);
    if exist(sid,'file') == 0
        iloop = iloop+1;
        continue
    else
        day_label{icount} = d(iloop).name;
        load(sid);
        load(stim);

        %re-organize into
        %day:stimulus_id:present_time_in_second:present_time_in_ms, row = day, col = trial
        for numcell = 1:length(stimulus_ids)
            response{icount,numcell} = horzcat(stimulus_ids{1,numcell}', stimulus_times{1,numcell}, stimulus_times{1,numcell} * 1000);
            sorted_response{icount,numcell} = sortrows(response{icount,numcell},1,'ascend'); %sort according to stimulus index
        end
        fprintf(num2str(iloop));
        iloop = iloop+1;
        icount = icount+1;
    end
end

save(fullfile(output_path, out_name), 'sorted_response','day_label')