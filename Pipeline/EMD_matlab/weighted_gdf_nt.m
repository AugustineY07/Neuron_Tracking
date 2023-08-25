function [E,L2] = weighted_gdf_nt(V1, V2, mw1, mw2, chan_pos, dim_mask, l2_weight)
%
% GDF   Ground distance between two vectors
%    [E] = GDF(F1, F2) is the ground distance between two feature vectors.
%
%    Example:
%    -------
%        v1 = [100, 40, 22];
%        v2 = [50, 100, 80];
%        ...
%        [e] = gdf(v1, v2);
%        ...
%
%    This file and its content belong to Ulas Yilmaz.
%    You are welcome to use it for non-commercial purposes, such as
%    student projects, research and personal interest. However,
%    you are not allowed to use it for commercial purposes, without
%    an explicit written and signed license agreement with Ulas Yilmaz.
%    Berlin University of Technology, Germany 2006.
%    http://www.cv.tu-berlin.de/~ulas/RaRF
%
% Modified from Ulas original by adding weights for each dimension and a 
% mask for which dimensions are included. 
% V1, V2 = [ndim,1], doubles, values for the two vectors

% hard coding these in the code for now, to be able to use original emd
% code and only change the distance function.
% dim_mask = [ndim,1], logical, to include or not
% dim_weights = weight for each, generally the standard deviation

% further modified to allow calling l2 distance between 2D waveforms from
% this function. 

dim_weights(1) = 0.1; % centroid x, 1/um^2
dim_weights(2) = 0.1; % centroid z, 1/um^2
dim_weights(3) = 1; % y from d1d2 data, 1/um^2
dim_weights(4) = 0.0025; % peak to peak amplitude, 1/uV^2
dim_weights(5) = 1000; % duration = peak-to-trough time, 1/msec^2
dim_weights(6) = 20000; % fwhm, 1/msec^2
dim_weights(7) = 7000; % peak to trough amplitude ratio, unitless
dim_weights(8) = 0.25; % pre peak amplitude, 1/uV^2
dim_weights(9) = 0.02; % vertical footprint 1/um^2
dim_weights(10) = l2_weight; % 'L2 norm' (estimated, need to run some bootstrapping)

% V1 = F1(i, 1:a);
% V2 = F2(j, 1:a);
% mw1 = squeeze(mw1(i,:,:));
% mw2 = squeeze(mw2(j,:,:));

V1 = V1(dim_mask(1:9));
V2 = V2(dim_mask(1:9));

w = dim_weights(dim_mask(1:9));
diffsq = (V2 - V1).^2;
if dim_mask(10) ~= 0
    % calculate l2 distance between these waveforms
    [~,wave_l2] = calcCenteredCorrL2(mw1,mw2, chan_pos, 5);
    diffsq = [diffsq,wave_l2^2];
    w = [w,dim_weights(10)]; %weights for all parameters
end

E = sqrt(dot(diffsq,w));
[~,L2] = calcCenteredCorrL2(mw1,mw2, chan_pos, 5);

end