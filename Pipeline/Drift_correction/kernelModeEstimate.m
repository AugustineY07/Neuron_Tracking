function modeEst = kernelModeEstimate(diffZ,varargin)

% Function to find the approximate mode given data. Based on the kernel
% approximation of the data density. 
%
% 
% 2022 - Adam Charles

data = diffZ;
if isvector(data)
    [fden, xden] = ksdensity(data(:)); % Compute one or two dimensional kernel density or distribution estimate, fden is the vector of density values. xden is the set of 100 (or 900) points.
    modeEst = xden(find(fden==max(fden),1,'first'));
elseif ismatrix(data)
    modeEst = zeros(1,size(data,2));
    parfor ll = 1:size(data,2)
        [fden, xden] = ksdensity(data(:,ll));
        modeEst(ll) = xden(find(fden==max(fden),1,'first'));
    end
elseif ndims(data) == 3
    modeEst = zeros(size(data,2),size(data,2));
    for ll = 1:size(data,1)
        parfor kk = 1:size(data,2)
            [fden, xden] = ksdensity(squeeze(data(ll,kk,:)));
            modeEst(ll,kk) = xden(find(fden==max(fden),1,'first'));
        end
        fprintf('.')
    end
    fprintf('\n')
end

end