function ks_call = readKS2label(fullpath, nUnit)

fid = fopen(fullpath,'r');
% allocate excess space to hold calls
% Depending on format, the number of calls may be < nUnit
ks_call = zeros([nUnit,1], 'logical');
% read line (this will be the header)
tline = fgetl(fid);
% read another line (first unit entry)
tline = fgetl(fid);
while ischar(tline)
    ln = split(tline);
    label = str2num(ln{1});
    call_str = ln{2};
    if strcmp(call_str,'good')
        ks_call(label+1) = 1;
    end
    tline = fgetl(fid);
end

ks_call = ks_call(1:nUnit);
%fprintf('%d out of %d units called good\n', sum(ks_call), nUnit);
end