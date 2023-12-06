function [f, L2] = gdm_nt(F1, F2,  mw1, mw2, chan_pos, dim_mask, l2_weight, xStep, zStep, Func)
% NT version modified from original to pass along mean waveform
% and channel position information to the ground distance function
% While technically just part of the signature, the need for the chan_pos
% information makes this the simpler format. Func still returns a single
% number.
%
% GDM   Ground distance matrix between two signatures
%    [F] = GDM(F1, F2, FUNC) is the ground distance matrix between
%    two signatures whose feature vectors are given in F1 and F2.
%    FUNC is a function which computes the ground distance between
%    two feature vectors.
%
%    Example:
%    -------
%        f1 = [[100, 40, 22]; [211, 20, 2]; [32, 190, 150]; [2, 100, 100]];
%        f2 = [[0, 0, 0]; [50, 100, 80]; [255, 255, 255]];
%        ...
%        [f] = gdm(f1, f2, @gmf);
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

% F1 = f1;
% F2 = f2;
% dim_mask = dim_mask_physical;

% number and length of feature vectors
[m a] = size(F1);
[n a] = size(F2);

% ground distance matrix
for i = 1:m
    for j = 1:n
%         f(i, j) = weighted_gdf_nt(F1(i, 1:a), F2(j, 1:a), squeeze(mw1(i,:,:)), squeeze(mw2(j,:,:)), chan_pos, dim_mask, l2_weight);
%         [f(i, j) L2(i,j)] = Func(F1(i, 1:a), F2(j, 1:a), squeeze(mw1(i,:,:)), squeeze(mw2(j,:,:)), chan_pos, dim_mask, l2_weight);
%         fprintf('F1(%d) = %d, F2(%d) = %d\n', i, F1, j, F2);
        [f(i, j) L2(i,j)] = weighted_gdf_nt(F1(i, 1:a), F2(j, 1:a), squeeze(mw1(i,:,:)), squeeze(mw2(j,:,:)), chan_pos, dim_mask, l2_weight, xStep, zStep);

    end
end

% gdm in column-vector form
f = f';
f = f(:);

L2 = L2';
L2 = L2(:);

end