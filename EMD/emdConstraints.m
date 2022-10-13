function Y = emdConstraints(X,N,M,opt)

% Function to compute EMD constraints for use in TFOCS A(F) = b, 
% where F is the flow matrix and b = [z1, z2, min(sum(z1), sum(z2)] are the constraints.
%
% As per TFOCS requirements, depending on 'opt', this function returns
%   opt == 0 :: The "size" of A [# constraints, # flows]
%   opt == 1 :: X is treated as F and A(F) is computed
%   opt == 2 :: X is treated as b and F = A^T(b) is computed
%
% 2022 - Adam Charles

% b = constraints
% X = P = flow matrix
% Y = opt1,2,3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if opt == 0 
    Y = [N+M+1, N*M]; % N = rows of S1, F(N*M): day 1 --> 2, A(N+M+1): F --> constraints
elseif opt == 1
    X = reshape(X,[N,M]); % X = F, Y = A(F) = b, set the dimension of flow matrix
    Y = [vec(sum(X,1)); vec(sum(X,2)); sum(sum(X))]; % A(F) = b, constraint for every flow
elseif opt == 2 % X = b, Y = A^T(b)
    Y = vec((vec(X(1:N))*ones(1,M) + ones(N,1)*(vec(X((N+1):(N+M)))') + X(N+M+1))'); % F = A^T(b)
else
    error('Can only compute A and A^T.\n')
end

end


function z = vec(z)
z = z(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
