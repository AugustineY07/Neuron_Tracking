function meas_out = wave_metrics(mw, chan_pos, stage,z_mode)

[nUnit, nChan, nt] = size(mw);
pp_all = squeeze(max(mw,[],3)-min(mw,[],3));
[pp_unit, pk_chan] = max(pp_all,[],2);
% get background of pp -- this is always nonzero, since it includes the
% noise.
ppvalues = reshape(pp_all,[nUnit*nChan,1]);
edges = 0:2:1000;
ppDist = histcounts(ppvalues,edges);
[~,maxbin] = max(ppDist);
backVal = (edges(maxbin+1)+edges(maxbin))/2;
pp_all = pp_all - 2*backVal;
neg_val = (pp_all<0);
pp_all(neg_val) = 0;

pk_wave = zeros(nUnit,nt);
for i = 1:nUnit
    pk_wave(i,:) = squeeze(mw(i,pk_chan(i),:));
end
min_unit = min(pk_wave,[],2);

meas_out = zeros(nUnit,2);
meas_out(:,1) = pp_unit;
meas_out(:,2) = min_unit;

pk_wave_upsamp = resample(pk_wave,200,82,'Dimension',2);
time_conv = (200/82)*30000/1000; % to convert from points to ms
[~,nt_up] = size(pk_wave_upsamp);
timestamps = (1:nt_up)/time_conv;
window = 20;  % for slope measurement

for i = 1:nUnit
    cw = squeeze(pk_wave_upsamp(i,:));
    npts_up = numel(cw);
    [~,trough_idx] = min(cw);
    [~,peak_idx] = max(cw);
    bPos = (cw(peak_idx) > abs(cw(trough_idx)));
    back_level = mean(cw(1:5));

    % duration/peak to trough time
    if bPos
        % this is a postive going peak, search from peak index to the next
        % trough
        [~, loc_trough_idx] = min(cw(peak_idx:end));
        duration = timestamps(peak_idx + loc_trough_idx-1) - timestamps(peak_idx);
    else
        % search forward from trough to find peak
        [~,loc_peak_idx] = max(cw(trough_idx:end));
        duration = timestamps(trough_idx + loc_peak_idx-1) - timestamps(trough_idx);
    end
    
    % fwhm
    fwhm = 0;
    if bPos
        threshold = cw(peak_idx) * 0.5;
        thresh_crossing_1 = min(timestamps(cw(1:peak_idx)>threshold));
        tc2_ind = peak_idx + min(find(cw(peak_idx:end)<threshold))-1;
        thresh_crossing_2 = timestamps(tc2_ind);     
    else
        threshold = cw(trough_idx) * 0.5;
        thresh_crossing_1 = min(timestamps(cw(1:trough_idx)<threshold));
        tc2_ind = trough_idx + min(find(cw(trough_idx:end)>threshold))-1;
        thresh_crossing_2 = timestamps(tc2_ind);  
    end
    if ~isempty(thresh_crossing_2) && ~isempty(thresh_crossing_1)
        fwhm = thresh_crossing_2 - thresh_crossing_1;
    end

    % PT ratio
    PT_ratio = 0;
    if cw(trough_idx) ~= 0
        PT_ratio = abs(cw(peak_idx)/cw(trough_idx));
    end

    % height of pre-peak, for negative going spikes
    pre_peak = 0;
    if ~bPos
        pre_peak = max(cw(1:trough_idx)) - back_level;
    end

    % recovery slope
    recovery_slope = 0;
    if bPos
       if trough_idx > npts_up-20
           window = npts_up - trough_idx;
       end
       x = timestamps(trough_idx:trough_idx+window-1); 
       X = [x', ones(window,1)];
       lreg = regress(cw(peak_idx:peak_idx+window-1)',X);
       lreg(1) = -1*lreg(1);
    else   
        % fit recover after repolarization (peak down to baseline)
        if peak_idx > npts_up-20
           window = npts_up - peak_idx;
       end
        x = timestamps(peak_idx:peak_idx+window-1);
        X = [x', ones(window,1)];
        lreg = regress(cw(peak_idx:peak_idx+window-1)',X);
    end
    recovery_slope = lreg(1) * 1e-3; % convert to V/s
    


    
    % 2D waveform metrics
    
    pp_unit = squeeze(pp_all(i,:))';
            % centroid calculation
%             distsq = (chan_pos(:,1)-chan_pos(pk_chan(i),1)).^2 + (chan_pos(:,2)-chan_pos(pk_chan(i),2)).^2;
%             sites_to_use = (distsq < (36)^2);
%             norm = sum(pp_unit(sites_to_use));
%             if norm > 0
%                 centX = sum(chan_pos(sites_to_use,1).*pp_unit(sites_to_use))/norm;
%                 centY = sum(chan_pos(sites_to_use,2).*pp_unit(sites_to_use))/norm;
%             else
%                 centX = min(chan_pos(:,1));
%                 centY = min(chan_pos(:,2));
%             end

% Julien Boussard - style fit of peak-to-peak voltage vs position
% if background sub pp_unit > 60 uV, attempt a fit of the
% background subtracted pp_all
if max(squeeze(pp_all(i,:))) > 60
    fitvals = fit_loc(i, pp_all, chan_pos);
    fitX = fitvals(1);
    fitZ = fitvals(2);
    fitY = fitvals(3);
else
    fitX = chan_pos(pk_chan(i),1);
    fitZ = chan_pos(pk_chan(i),2);
    fitY = -1;      % a marker for no fit
end

    sp_thresh = 0.2*max(pp_unit);
    chan_above = pp_unit > sp_thresh;
    zmax = max(chan_pos(chan_above,2));
    zmin = min(chan_pos(chan_above,2));
    z_sp = 0;
    if ~isempty(zmax) && ~isempty(zmin)
        z_sp = zmax(1)-zmin(1);
    end
   

    % fill in meas_out for this unit
    meas_out(i,3) = duration;
    meas_out(i,4) = fwhm;
    meas_out(i,5) = PT_ratio;
    meas_out(i,6) = pre_peak;
    meas_out(i,7) = recovery_slope;
    meas_out(i,8) = z_sp;
    meas_out(i,9) = fitX;
    switch stage
        case 'pre'
            meas_out(i,10) = fitZ;
        case 'post'
            meas_out(i,10) = fitZ - z_mode;
    end
    meas_out(i,11) = fitY;

end

end
