function ks2_working_toBinary( varargin )

% Read data preprocessed by KS2 (whitenend, post datashift) and write out
% as a standard binary.

% IMPORTANT NOTE: the preprocessed file only contains channels that KS2
% uses for sorting! Channels with connected=0 in the chanMap file OR
% eliminated due to low spike count (ops.minfr_goodchannels > 0) will not
% be included. In that case, use either the info in the rez file or the phy
% output to specify the correct site geometry and channel map for sorting.

% 2nd important note:
% Make sure you set the correct version of kilosort. Supported versions are
% '2.0' and '2.5' and '3.0'
KSver = '2.5'
triggerStr = 'tcat'; % set to t0 for data not processed by CatGT, tcat for processed data
useRezFroc = 0; % set to 1 if the path to fproc in rez.mat is correct, 0 if it needs to be selected separately

lenTS = numel(triggerStr);

% make sure the rez file is pointing to the current preprocessed data file.
if isempty(varargin )
    [rezName, rezPath] = uigetfile('*.mat','Select rez file from KS sort');
    rezFullPath = fullfile( rezPath, rezName );
    modStr = 'new';
    [mmName, mmPath] = uigetfile('*.meta','Select original metadata file');
    modelMetaFullPath = fullfile(mmPath, mmName);
    if useRezFroc ~= 1     
        [fprocName, fprocPath] = uigetfile({'*.bin','*.dat'}, 'Select fproc file');
        fprocFullPath = fullfile(fprocPath,fprocName);
    end
        
else
    % called with rezFullPath
    inputCell = varargin(1);
    rezFullPath = inputCell{1};
    inputCell = varargin(2);
    modStr = inputCell{1};
    inputCell = varargin(3);
    modelMetaFullPath = inputCell{1};
    if useRezFroc ~= 1  
        inputCell = varargin(4);
        fprocFullPath = inputCell{1};
    end

end

% build output name, get copy of metaName
load( rezFullPath );
ops = rez.ops;

if useRezFroc
    fprocFullPath = ops.fproc;
end

[binPath, binName, binExt] = fileparts(rez.ops.fbinary);
tPos = strfind(binName, '.imec');
temp = extractBetween(binName,1,tPos-(lenTS+1));
baseName = temp{1};

suffix = extractAfter(binName,triggerStr);
outName = sprintf( '%s%s%s%s', baseName, modStr, suffix, binExt);
if useRezFroc
    % 'Typically' want the output to go in the directory with the input
    [outPath, ~, ~] = fileparts(binPath);
    outFullPath = fullfile(binPath, outName);
else
    % 'Typically' working from a copy of the processed file, and the 
    % unwhitened file goes there
    [outPath, ~, ~] = fileparts(fprocFullPath);
    outFullPath = fullfile(outPath,outName);
end
metaName = sprintf( '%s%s', extractBefore(outName,'bin'), 'meta');

% KS2 creates a temporary whitened file consisting of batches of NT time
% points stored as NT rows by NChan columns (so each column the time trace
% for a single channel). These blocks are overlapped in time by nt.buff
% timepoints, which get trimmed off in analysis to avoid some filtering
% artifacts. To translate these data to "standard binary format, need to
% read each block, trim off the buffer, unwhiten and rescale, and finally
% tranpose to Nchan rows by NT-ntbuff columns.

% KS2.5 aims to make a readable binary from it's datashifted data. Rather
% than overlapping the batches, it reads in some extra points for filtering
% and then trims them back off. ops.ntbuff = ntb
% read NTbuff = NT + 3*ntb points
% for a standard batch (neither first nor last) read:
%       --ntb points overlapping the last batch
%       --NT points that belong to this batch
%       --2*ntb more points; first set of ntb will be blended with next
%       batch, 2nd ntb just filtering buffer
% After fitering, points ntb+1:ntb are blended with "dat_prev" which is
% NT+ntb+1:NT+2*ntb saved from the previous batch.
% Batch zero gets ntb zeros prepended to its data and blended with the
% initialized dat_prev (also zeros). 
% After filtering, the data is whitened and then transposed bacyk to NChan
% rows by NT columns to save. When these batches are read for sorting in
% learnTemplates, the data is transposed after reading (new in KS25)




fid = fopen(fprocFullPath, 'r');
fidW = fopen(outFullPath, 'w'); % open for writing processed data, transposed

if strcmp(KSver,'2.0')
    batchstart = 0:ops.NT:ops.NT*ops.Nbatch; % batches start at these timepoints
    for ibatch = 1:ops.Nbatch
        offset = 2 * ops.Nchan*batchstart(ibatch); % binary file offset in bytes
        fseek(fid, offset, 'bof');
        dat = fread(fid, [ops.NT ops.Nchan], '*int16');
        % Due to clumsy arithmetic in KS, the first two batches overlap by 
        % 2*ops.ntbuff, while later batches overlap by 1*ntbuff
        % trim samples of the end of the current batch accordingly
        if ibatch == 1
            dat = dat(1:ops.NT-2*ops.ntbuff, :);
        else
            dat = dat(1:ops.NT-ops.ntbuff, :);
        end
        % rescale by the average of the diagonal of the whitening
        % matrix. This will keep the voltage values in range during the 2nd
        % sort.
        %norm = mean(diag(rez.Wrot));
        %dat = int16(single(dat)./norm);
        dat = int16(single(dat)/rez.Wrot);
        fwrite(fidW, dat', 'int16'); % write transposed batch to binary, these chunks are in order
    end
elseif strcmp(KSver,'2.5') || strcmp(KSver,'3.0')
    %Batches already stored end to end and transposed. Just need to read,
    %transpose, unwhiten, tranpose back and store
    batchstart = 0:ops.NT:ops.NT*ops.Nbatch; % batches start at these timepoints
    for ibatch = 1:ops.Nbatch
        offset = 2 * ops.Nchan*batchstart(ibatch); % binary file offset in bytes
        fseek(fid, offset, 'bof');
        dat = fread(fid, [ops.Nchan ops.NT], '*int16');
        % if skipping unwhitening, rescale by the average of the diagonal of the whitening
        % matrix. This will keep the voltage values in range during the 2nd
        % sort.
        %norm = mean(diag(rez.Wrot));
        %dat = int16(single(dat)./norm);
        % if unwhitening, skip the rescale step. Transpose to NT by Nchan
        % before unwhitening
        dat = int16(single(dat')/rez.Wrot);
        fwrite(fidW, dat', 'int16'); % write transposed batch to binary, these chunks are in order
    end
else
    fprintf( 'Unknown version of Kilsort.\n');

end

fclose(fid);
fclose(fidW);

if strlength(modelMetaFullPath) > 0
    
    fp = dir(outFullPath);
    newTags = cell(4,1);
    newTag{1} = sprintf('%s%d', 'fileSizeBytes=', fp.bytes);
    newTag{2} = sprintf('%s%d', 'nSavedChans=', ops.Nchan);
    newTag{3} = sprintf('%s%d%s', 'snsApLfSy=', ops.Nchan, ',0,0');
    newTag{4} = sprintf('%s%d', 'snsSaveChanSubset=0:',ops.Nchan-1);
    repTags = cell(4,1);
    repTags{1} = 'fileSizeBytes';
    repTags{2} = 'nSavedChans';
    repTags{3} = 'snsApLfSy';
    repTags{4} = 'snsSaveChanSubset';
    
    fmodel = fopen( modelMetaFullPath, 'r');
    fmeta = fopen( fullfile(outPath, metaName), 'w');
    
    tline = fgetl(fmodel);
    while ischar(tline)
        currTag = extractBefore(tline,'=');
        tagFound = find(strcmp(repTags, currTag));
        if isempty(tagFound)
            %copy over this line as is
            fprintf(fmeta, '%s\n', tline );
        else
            fprintf('found: %s\n', repTags{tagFound} );
            fprintf(fmeta, '%s\n', newTag{tagFound} );
        end  
        tline = fgetl(fmodel);
    end
    fclose(fmeta);
    fclose(fmodel);
end
fprintf( 'Output file has %d channels\n', ops.Nchan );

end
