function ks_call = readKS2label(fullpath)

fid = fopen(fullpath,'r');
nUnit = 0;
%allocate excess space to hold ks calls
ks_call = zeros([5000,1], 'logical');

% read line (this will be the header)
tline = fgetl(fid);
% read another line (first unit entry)
tline = fgetl(fid);
while ischar(tline)
    nUnit = nUnit + 1;
    ln = split(tline);
    call_str = ln{2};
    if strcmp(call_str,'good')
        ks_call(nUnit) = 1;
    end
    tline = fgetl(fid);
end

ks_call = ks_call(1:nUnit);
fprintf('%d out of %d units called good\n', sum(ks_call), nUnit);
end