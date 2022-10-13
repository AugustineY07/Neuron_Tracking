function [P,cost] = runTfocsMatching(f1,f2,optSel)



mu = 1;
opts.nonnegativity = true;                                                 % Make sure F is non-negative
opts.tol = 1.0000e-08;
opts.printStopCrit = true;

switch lower(optSel)
    case {'constrained','both'}
        for iclu = 1:size(f1,1)
            for jclu = 1:size(f2,1)
                x = abs(f1(iclu,1)-f2(jclu,1));
                z = abs(f1(iclu,2)-f2(jclu,2));
                if x<=2 && z<=2
                    C(iclu,jclu) = sqrt((f1(iclu,1)-f2(jclu,1))^2 + (f1(iclu,2)-f2(jclu,2))^2 + (f1(iclu,3)-f2(jclu,3))^2);
                else
                    C(iclu,jclu) = 0;
                end
            end
        end
        C = C';
        C = C(:); % vectorized vertically
        ind = find(C <= inf);
        C_new = C(ind);
    otherwise
        C = gdm(f1, f2, @gdf);
end

%%
switch lower(optSel)
    case {'partial','both'}
        tau = 0.8; %partial parameter
    otherwise
        tau = 1;
end
b = [ones(size(f1,1),1); ones(size(f2,1),1); round(tau*min(size(f1,1), size(f2,1)))]; % parameter that A computed: all S1, all S2, smaller of total weights between S1 and S2

switch lower(optSel)
    case {'constrained','both'}
        A = @(X,T) emdConstraints_threshold(X, ind, size(f1,1),size(f2,1),T); % Constraints matrix, T is option, takes value 0, 1, and 2
        res  = solver_sLP(C_new, A, b, mu, [], [], opts); % Run the Generic linear programming in standard form. c and b are vectors, A is matrix
        cost = sum(res.*C_new);                                                        % EMD value
        P = zeros(size(f1,1),size(f2,1));
        P(ind) = res;
    otherwise
        C = gdm(f1, f2, @gdf);
        A = @(X,T) emdConstraints(X,size(f1,1),size(f2,1),T); % Constraints matrix, T is option, takes value 0, 1, and 2

        res = solver_sLP(C, A, b, mu, [], [], opts);                               % Run the Generic linear programming in standard form. c and b are vectors, A is matrix
        cost = sum(res.*C);                                                        % EMD value
        P = reshape(res,size(f2,1),size(f1,1));                                    % Reshape output
end

                       
end
