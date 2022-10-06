% This function contains all EMD calculation with all options

clear all;
close all;

pars.nTrials = 10;
pars.npts = 50; %number of points per trial
pars.xRange = [-50,50];
pars.zRange = [-50,770];
pars.yRange = [-10,10];
pars.y_drift = 20;
pars.pos_err_std = 5;
pars.x_err_std = 1;
pars.y_err_std = 1;
pars.z_err_std = 1;
pars.plot = 1;
pars.save = 1; %0 = not save, 1 = save data

pars.rootD = 'C:\Users\labadmin\Desktop\EMD_test_simulator\'; % the simulated data files are in this folder
pars.alg = 3; %1 = matlab, 2 = CVX, 3 = TFOCS
pars.runType = 1; %0 = data generation, 1 = run EMD
pars.dimension = 3; %2 = 2D, 3 = 3D
pars.errorType = 2; %1 = pos_err_std, 2 = x_y_z_error
pars.mode = 2; %0 = standard, 1 = lumpiness, 2 = gainloss
pars.block = 4; %lumpiness: number of blocks
pars.changeType = 1; %gainloss: 1 = gain, 2 = loss, 3 = mixed
pars.change = [1:1:5]; %number of lost/gain points



%
switch pars.runType %data v.s run
    case 0
        data_generator(pars);
    case 1
        numIncorrect = [];
        perIncorrect = [];
        switch pars.dimension
            case 2
                dim = '2D';
            case 3
                dim = '3D';
        end
        switch pars.errorType
            case 1
                errorType = 'Pos';
            case 2
                errorType = 'xyz';
        end
        switch pars.mode
            case 0
                mode = 'standard';
                name = ['C:\Users\labadmin\Desktop\EMD_test_simulator\data\', dim,'\', errorType,'\', mode,'\'];
            case 1
                mode = 'lumpiness';
                name = ['C:\Users\labadmin\Desktop\EMD_test_simulator\data\', dim,'\', errorType,'\', mode,'\'];
            case 2
                mode = 'gainloss';
                switch pars.changeType
                    case 1
                        changeType = 'gain';
                    case 2
                        changeType = 'loss';
                    case 3
                        changeType = 'mixed';
                end
                name = ['C:\Users\labadmin\Desktop\EMD_test_simulator\data\', dim,'\', errorType,'\', mode,'\', changeType,'\'];
        end

        for ich = 1:length(pars.change)
            for i = 1:pars.nTrials
                %points = load([name,'trial', num2str(i), ',', ' drift=', num2str(pars.y_drift), ' error=', num2str(pars.x_err_std), '&', num2str(pars.y_err_std), '&', num2str(pars.z_err_std), '.mat']);
                points = load([pars.rootD, mode, '\data\', changeType, '\trial', num2str(i), ',', 'npts', num2str(pars.npts), changeType, num2str(pars.change(ich)), ',error=', num2str(pars.x_err_std), '&', num2str(pars.y_err_std),'&', num2str(pars.z_err_std), '.mat']);

                f1 = points.f1;
                f2 = points.f2;
                f2_more = points.f2_more;
                loss = points.loss;

                % set weights = 1
                w1 = ones(length(f1),1);%ones([length(f1),1])/length(f1);
                w2 = ones(length(f2),1);%ones([length(f2),1])/length(f2);

                switch pars.alg
                    case 1 %matlab
                        mFolder = 'matlab';
                        [x, fval] = emd(f1, f2, w1, w2, @gdf);
                        P = reshape(x,[length(f2),length(f1)]);
                        C = gdm(f1, f2, @gdf); % Distance matrix
                    case 2 %CVX
                        mFolder = 'CVX';
                        N0 = length(f1);
                        N = length(f2);
                        C = [];
                        %b = min(min(f1),min(f2)); %shift for any negative values, determined by the smallest point
                        for iclu = 1:N
                            for jclu = 1:N0
                                C(iclu,jclu) = sqrt((f2(iclu,1)-f1(jclu,1))^2 + (f2(iclu,2)-f1(jclu,2))^2 + (f2(iclu,3)-f1(jclu,3))^2);
                            end
                        end
                        % if b < 0
                        %     x = x + abs(b);
                        %     x = (x - min(x)) / max(x - min(x)); %normalize x
                        %     x0 = x0 + abs(b);
                        %     x0 = (x0 - min(x0)) / max(x0 - min(x0));
                        % end

                        % Solve
                        tic
                        cvx_begin quiet
                        % cvx_solver
                        cvx_solver SDPT3 %optimization solver, if slow, try smaller x first
                        cvx_precision high
                        variable P(N,N0) nonnegative; %f_ij >= 0
                        minimize(C(:)'*P(:));
                        subject to
                        sum(P,2)  <= w2(:); %can't move more to destination
                        sum(P,1)' <= w1(:); %can't move more from origin
                        sum(P(:)) == min(sum(w2),sum(w1)); %total amount is the smaller mass
                        cvx_end

                        stats.emd     =   sum(C(:)'*P(:));
                        stats.runtime = toc;
                        stats.nbr_iter = cvx_slvitr; % the number of iterations taken by the solver

                        EMD_cost(ich,i) = stats.emd;
                        Pmatrix{ich,i} = P;

                    case 3 %TFOCS
                        mFolder = 'TFOCS';
                        C = [];
                        mu = 1;
                        for iclu = 1:size(f1,1)
                            for jclu = 1:size(f2,1)
                                C(iclu,jclu) = sqrt((f1(iclu,1)-f2(jclu,1))^2 + (f1(iclu,2)-f2(jclu,2))^2 + (f1(iclu,3)-f2(jclu,3))^2);
                            end
                        end
                        C = C(:); % vectorized vertically
                        A = @(X,T) emdConstraints(X,size(f1,1),size(f2,1),T); % Constraints matrix, T is option, takes value 0, 1, and 2
                        b = [ones(size(f1,1),1); ones(size(f2,1),1); min(size(f1,1), size(f2,1))]; % parameter that A computed: all S1, all S2, smaller of total weights between S1 and S2

                        opts.nonnegativity = true;                                                 % Make sure F is non-negative
                        res = solver_sLP(C, A, b, mu, [], [], opts);                               % Run the Generic linear programming in standard form. c and b are vectors, A is matrix
                        cost = sum(res.*C);                                                        % EMD value
                        P = reshape(res,size(f1,1),size(f2,1));                                    % Reshape output

%                         mu = 1;
%                         for iclu = 1:size(f1,1)
%                             for jclu = 1:size(f2,1)
%                                 x = abs(f1(iclu,1)-f2(jclu,1));
%                                 z = abs(f1(iclu,2)-f2(jclu,2));
%             %                     if x<=1 && z<=1
%                                     C(iclu,jclu) = sqrt((f1(iclu,1)-f2(jclu,1))^2 + (f1(iclu,2)-f2(jclu,2))^2 + (f1(iclu,3)-f2(jclu,3))^2);
%             %                     else
%             %                         C(iclu,jclu) = 0;
%             %                     end
%                             end
%                         end
%                         C = C';
%                         C = C(:); % vectorized vertically
%                         tau = 1; %partial parameter
%                         b = [ones(size(f1,1),1); ones(size(f2,1),1); round(tau*min(size(f1,1), size(f2,1)))]; % parameter that A computed: all S1, all S2, smaller of total weights between S1 and S2
%                         opts.nonnegativity = true;                                                 % Make sure F is non-negative
%                         opts.printStopCrit = true;
%                         %opts.restart =
%                         %opts.maxCounts = [ Inf, Inf, iteration, Inf, Inf ];
%                         ind = find(C <= 20);
%                         C_new = C(ind);
%                         A = @(X,T) emdConstraints_threshold(X, ind, size(f1,1),size(f2,1),T); % Constraints matrix, T is option, takes value 0, 1, and 2
%                         res  = solver_sLP_Rplus(C_new, A, b, mu, [], [], opts); % Run the Generic linear programming in standard form. c and b are vectors, A is matrix
            
%                         cost = sum(res.*C_new);                                                        % EMD value
%                         P = zeros(size(f2,1),size(f1,1));
%                         P(ind) = res;



                        EMD_cost(ich,i) = cost;
                        Pmatrix{ich,i} = P;
                end
                % calculate norm
                C = reshape(C,[length(f1),length(f2)]);
                %distance(ip,i) = mean(C,'all'); %calculate mean interneuron distance
                sortC = sort(C);
                neighborAve(ich,i) = mean(sortC(2,:)); %find mean of all nearest neighbor distance
                move = diag(C);
                offsetMean(ich,i) = mean(move); %calculate mean offset between points

                noMatch = [];
                switch pars.changeType
                    case 1 %gain
%                         % pad added points with 0 in P cols
%                         ind_change = ones(1,length(f2));
%                         ind_change(length(f1)+1:length(f1)+ich) = 0;
%                         ind = find(ind_change);
%                         P_new = zeros(length(f2), length(f2));
%                         P_new(:,ind) = P;
% 
%                         % for each row, find the nonzero max column
%                         P_new(P_new==0) = nan;
%                         [~,max_ind] = max(P_new,[],1);
% 
%                         % match index in f1
%                         f1_full = zeros(length(f2),3);
%                         f1_full(ind,:) = f1; %recover f1 after change
% 
%                         % reorganize matches
%                         f2_match = f2(max_ind,:); %matched points in order
%                         numIncorrect(ich, i) = sum(max_ind ~= [1:sum(length(f2))]);
%                         perIncorrect(ich, i) = sum(max_ind ~= [1:sum(length(f2))])/sum(length(f2));
% 
%                         % calculate distance of matched pairs
%                         match_dis = zeros(length(f2),2);
%                         match_dis(:,1) = max_ind;
%                         count = 0;
%                         for ir = 1:length(f2)
%                             if ~ismember(ir, f2_more)
%                                 count = count + 1;
%                                 match_dis(count,2) = sqrt( (f1_full(ir,1)-f2(match_dis(ir),1))^2 + (f1_full(ir,2)-f2(match_dis(ir),2))^2 + (f1_full(ir,3)-f2(match_dis(ir),3))^2 );
%                                 correct_dis(count,1) = sqrt((f1_full(ir,1)-f2(ir,1))^2 + (f1_full(ir,2)-f2(ir,2))^2 + (f1_full(ir,3)-f2(ir,3))^2);
%                             end
%                         end
%                         figure(i)
%                         histogram(match_dis(:,2),10)
%                         figure(10+i)
%                         histogram(correct_dis)




                        P_new = zeros(length(f1), length(f2));
                        P_new = P;
                        P_new(P_new==0) = nan;
                        noMatch(:,1) = find(all(isnan(P_new),2)); %find rows with all nan
                        [~,max_ind] = max(P_new,[],2); %nonzero max in every row
                        % match index
                        f2_same = points.f2_same; %points preserved in f2
                        [q, f2_unchanged] = ismember(f2_same, f2, 'rows'); %the new index in f2 corresponds to true f1 pair index

                        % reorganize matches
                        f2_match = f2(max_ind,:); %f2 index 
                        numIncorrect(ich, i) = sum(max_ind(:) ~= f2_unchanged);
                        perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1)+ich);
                        if ~isempty(noMatch) %number of correct unmatch
                            correct_no = sum(ismember(noMatch,loss)) + sum(ismember([length(f1)+1:1:length(f1)+ich],noMatch));
                            numIncorrect(ich, i) = numIncorrect(ich, i) - correct_no;
                            perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1)+ich);
                        end

                        % calculate distance of matched pairs
                        match_dis = zeros(length(f1),4);
                        match_dis(:,1) = 1:length(f1); %f1 index
                        match_dis(:,2) = f2_unchanged(1:length(f1)); %f2 index corresponds to f1
                        match_dis(:,3) = max_ind(1:length(f1)); %result match index in f2 scale
                        correct_dis = zeros(length(f1),1);
                        for ir = 1:length(f1)
                            match_dis(ir,4) = sqrt( (f1(ir,1)-f2(match_dis(ir,3),1))^2 + (f1(ir,2)-f2(match_dis(ir,3),2))^2 + (f1(ir,3)-f2(match_dis(ir,3),3))^2 );
                            correct_dis(ir) = sqrt( (f1(ir,1)-f2(match_dis(ir,2),1))^2 + (f1(ir,2)-f2(match_dis(ir,2),2))^2 + (f1(ir,3)-f2(match_dis(ir,2),3))^2 );
                        end
                        figure(i)
                        histogram(match_dis(:,4),10)
                        figure(10+i)
                        histogram(correct_dis)



                    case 2 %loss
%                         % pad missing points with 0 in P rows
%                         ind_change = ones(1,length(f1));
%                         ind_change(loss) = 0;
%                         ind = find(ind_change);
%                         P_new = zeros(length(f1), length(f1));
%                         P_new(ind,:) = P;
% 
%                         % for each row, find the nonzero max column
%                         P_new(P_new==0) = nan;
%                         [~,max_ind] = max(P_new,[],1);
% 
%                         % match index in f2
%                         f2_full = zeros(length(f1),3);
%                         f2_full(ind,:) = f2; %recover f2 before change
% 
%                         % reorganize matches
%                         f2_match = f2_full(max_ind,:); %matched points in order
%                         numIncorrect(ich, i) = sum(max_ind ~= [1:sum(pars.npts)]);
%                         perIncorrect(ich, i) = sum(max_ind ~= [1:sum(pars.npts)])/sum(pars.npts);
% 
%                         % calculate distance of matched pairs
%                         match_dis = zeros(length(f1),2);
%                         match_dis(:,1) = max_ind;
%                         for ir = 1:length(f1)
%                             match_dis(ir,2) = sqrt( (f1(ir,1)-f2_full(match_dis(ir),1))^2 + (f1(ir,2)-f2_full(match_dis(ir),2))^2 + (f1(ir,3)-f2_full(match_dis(ir),3))^2 );
%                             correct_dis(ir,1) = sqrt((f1(ir,1)-f2_full(ir,1))^2 + (f1(ir,2)-f2_full(ir,2))^2 + (f1(ir,3)-f2_full(ir,3))^2);
%                         end
%                         figure(i)
%                         histogram(match_dis(:,2),10)
%                         figure(10+i)
%                         histogram(correct_dis)

                        P_new = zeros(length(f1), length(f2));
                        P_new = P;
                        P_new(P_new==0) = nan;
                        noMatch(:,1) = find(all(isnan(P_new),2)); %find rows with all nan
                        [~,max_ind] = max(P_new,[],2); %nonzero max in every row
                        % match index
                        f2_same = points.f2_same; 
                        [q, f2_unchanged] = ismember(f2_same, f2, 'rows'); %the new index in f2 corresponds to true f1 pair index

                        % reorganize matches
                        f2_match = f2(max_ind,:); %f2 index 
                        numIncorrect(ich, i) = sum(max_ind(:) ~= f2_unchanged);
                        perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1));
                        if ~isempty(noMatch) %number of correct unmatch
                            correct_no = sum(ismember(noMatch,loss)) + sum(ismember([length(f1)+1:1:length(f1)+ich],noMatch));
                            numIncorrect(ich, i) = numIncorrect(ich, i) - correct_no;
                            perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1));
                        end

                        % calculate distance of matched pairs
                        match_dis = zeros(length(f1),4);
                        match_dis(:,1) = 1:length(f1); %f1 index
                        match_dis(:,2) = f2_unchanged; %f2 index corresponds to f1
                        match_dis(:,3) = max_ind; %result match index in f2 scale
                        correct_dis = zeros(length(f1),1);
                        for ir = 1:length(f1)
                            match_dis(ir,4) = sqrt( (f1(ir,1)-f2(match_dis(ir,3),1))^2 + (f1(ir,2)-f2(match_dis(ir,3),2))^2 + (f1(ir,3)-f2(match_dis(ir,3),3))^2 );
                            if ~ismember(ir,loss)
                                correct_dis(ir) = sqrt( (f1(ir,1)-f2(match_dis(ir,2),1))^2 + (f1(ir,2)-f2(match_dis(ir,2),2))^2 + (f1(ir,3)-f2(match_dis(ir,2),3))^2 );
                            else
                                correct_dis(ir) = sqrt( (f1(ir,1)-0)^2 + (f1(ir,2)-0)^2 + (f1(ir,3)-0)^2 );
                            end
                        end
                        figure(i)
                        histogram(match_dis(:,4),10)
                        figure(10+i)
                        histogram(correct_dis)


                    case 3 %mixed
                        P_new = zeros(length(f1), length(f2));
                        P_new = P;
                        P_new(P_new==0) = nan;
                        noMatch(:,1) = find(all(isnan(P_new),2)); %find rows with all nan
                        [~,max_ind] = max(P_new,[],2); %nonzero max
                        % match index
                        f2_same = points.f2_same;
                        [q, f2_unchanged] = ismember(f2_same, f2, 'rows'); %the new index in f2 preserved from f1 

                        % reorganize matches
                        f2_match = f2(max_ind,:); %f2 index 
                        numIncorrect(ich, i) = sum(max_ind(1:length(f1)) ~= f2_unchanged');
                        perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1)+ich);
                        if ~isempty(noMatch) %number of correct unmatch
                            correct_no = sum(ismember(noMatch,loss)) + sum(ismember([length(f1)+1:1:length(f1)+ich],noMatch));
                            numIncorrect(ich, i) = numIncorrect(ich, i) - correct_no;
                            perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1)+ich);
                        end

                        % calculate distance of matched pairs
                        match_dis = zeros(length(f1),4);
                        match_dis(:,1) = 1:length(f1); %f1 index
                        match_dis(:,2) = f2_unchanged; %f2 index corresponds to f1
                        match_dis(:,3) = max_ind; %result match index in f2 scale
                        correct_dis = zeros(length(f1),1);
                        for ir = 1:length(f1)
                            match_dis(ir,4) = sqrt( (f1(ir,1)-f2(match_dis(ir,3),1))^2 + (f1(ir,2)-f2(match_dis(ir,3),2))^2 + (f1(ir,3)-f2(match_dis(ir,3),3))^2 );
                            if ~ismember(ir,loss)
                                correct_dis(ir) = sqrt( (f1(ir,1)-f2(match_dis(ir,2),1))^2 + (f1(ir,2)-f2(match_dis(ir,2),2))^2 + (f1(ir,3)-f2(match_dis(ir,2),3))^2 );
                            else
                                correct_dis(ir) = sqrt( (f1(ir,1)-0)^2 + (f1(ir,2)-0)^2 + (f1(ir,3)-0)^2 );
                            end
                        end
                        figure(i)
                        histogram(match_dis(:,4),10)
                        figure(10+i)
                        histogram(correct_dis)
                end
            end



        end
end

% numIncorrect_cvx_gain = numIncorrect;
% save 'C:\Users\labadmin\Desktop\EMD_test_simulator\Algorithm comparison\gain\cvx_gain50.mat' numIncorrect_cvx_gain


% ----------------------- Helper functions ----------------------------------------
% Data generator
% function data_generator(pars)
%
% end