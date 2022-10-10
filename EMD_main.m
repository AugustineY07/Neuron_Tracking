% This function contains all EMD calculation with all options

clear all;
close all;

addpath(genpath('C:\Users\labadmin\Desktop\EMD_test_simulator')) % path to simulator folder
addpath(genpath('D:\Data\Pipeline\EMD\Localization_EMD_TFOCS')) % TFOCS package

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
pars.alg = 'TFOCS'; %1 = matlab, 2 = CVX, 3 = TFOCS
pars.runType = 1; %0 = data generation, 1 = run EMD
pars.dimension = '3D'; %2 = '2D', 3 = '3D'
pars.errorType = 'xyz'; %'Pos' = pos_err_std, 'xyz' = x_y_z_error
pars.mode = 'gainloss'; %0 = standard, 1 = lumpiness, 2 = gainloss
pars.block = 4; %lumpiness: number of blocks
pars.changeType = 'mixed'; %gainloss: 1 = gain, 2 = loss, 3 = mixed
pars.change = [1:1:5]; %number of lost/gain points



% switch pars.dimension
%     case 2
%         dim = '2D';
%     case 3
%         dim = '3D';
% end
% switch pars.errorType
%     case 1
%         errorType = 'Pos';
%     case 2
%         errorType = 'xyz';
% end
switch pars.mode
    case {'standard','lumpiness'}
        name = ['C:\Users\labadmin\Desktop\EMD_test_simulator\data\', pars.dimension,'\', pars.errorType,'\', pars.mode,'\'];
    case 'gainloss'
        name = ['C:\Users\labadmin\Desktop\EMD_test_simulator\data\', pars.dimension,'\', pars.errorType,'\', pars.mode,'\', pars.changeType,'\'];
end

switch pars.runType %data v.s run
    case 0
        data_generator(pars, name);
    case 1
        numIncorrect = [];
        perIncorrect = [];

        for ich = 1:length(pars.change)
            for i = 1:pars.nTrials
                %points = load([name,'trial', num2str(i), ',', ' drift=', num2str(pars.y_drift), ' error=', num2str(pars.x_err_std), '&', num2str(pars.y_err_std), '&', num2str(pars.z_err_std), '.mat']);
                points = load([pars.rootD, pars.mode, '\data\', pars.changeType, '\trial', num2str(i), ',', 'npts', num2str(pars.npts), pars.changeType, num2str(pars.change(ich)), ',error=', num2str(pars.x_err_std), '&', num2str(pars.y_err_std),'&', num2str(pars.z_err_std), '.mat']);

                f1 = points.f1;
                f2 = points.f2;
                f2_more = points.f2_more;
                loss = points.loss;

                % set weights = 1
                w1 = ones(length(f1),1);%ones([length(f1),1])/length(f1);
                w2 = ones(length(f2),1);%ones([length(f2),1])/length(f2);

                C = [];

                switch pars.alg
                    case 'matlab'
                        [x, fval] = emd(f1, f2, w1, w2, @gdf);
                        P = reshape(x,[length(f2),length(f1)]);
                        C = gdm(f1, f2, @gdf); % Distance matrix

                        EMD_cost(ich,i) = fval;
                        Pmatrix{ich,i} = P;
                    case 'CVX'
                        %b = min(min(f1),min(f2)); %shift for any negative values, determined by the smallest point
                        for iclu = 1:length(f1)
                            for jclu = 1:length(f2)
                                C(iclu,jclu) = sqrt((f1(iclu,1)-f2(jclu,1))^2 + (f1(iclu,2)-f2(jclu,2))^2 + (f1(iclu,3)-f2(jclu,3))^2);
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
                        variable P(length(f1),length(f2)) nonnegative; %f_ij >= 0
                        minimize(C(:)'*P(:)); %vectorize C by col, C became f2*f1
                        subject to
                        sum(P,2)  <= w1(:); %can't move more to destination
                        sum(P,1)' <= w2(:); %can't move more from origin
                        sum(P(:)) == min(sum(w2),sum(w1)); %total amount is the smaller mass
                        cvx_end

                        stats.emd     =   sum(C(:)'*P(:));
                        stats.runtime = toc;
                        stats.nbr_iter = cvx_slvitr; % the number of iterations taken by the solver

                        EMD_cost(ich,i) = stats.emd;
                        Pmatrix{ich,i} = P;

                    case 'TFOCS'
                        %                         mu = 1;
                        %                         for iclu = 1:size(f1,1)
                        %                             for jclu = 1:size(f2,1)
                        %                                 C(iclu,jclu) = sqrt((f1(iclu,1)-f2(jclu,1))^2 + (f1(iclu,2)-f2(jclu,2))^2 + (f1(iclu,3)-f2(jclu,3))^2);
                        %                             end
                        %                         end
                        %                         C = C(:); % vectorized by col
                        %                         A = @(X,T) emdConstraints(X,size(f1,1),size(f2,1),T); % Constraints matrix, T is option, takes value 0, 1, and 2
                        %                         b = [ones(size(f1,1),1); ones(size(f2,1),1); min(size(f1,1), size(f2,1))]; % parameter that A computed: all S1, all S2, smaller of total weights between S1 and S2
                        %
                        %                         opts.nonnegativity = true;                                                 % Make sure F is non-negative
                        %                         res = solver_sLP(C, A, b, mu, [], [], opts);                               % Run the Generic linear programming in standard form. c and b are vectors, A is matrix
                        %                         cost = sum(res.*C);                                                        % EMD value
                        %                         P = reshape(res,size(f1,1),size(f2,1));                                    % Reshape output

                        mu = 1;
                        for iclu = 1:size(f1,1)
                            for jclu = 1:size(f2,1)
                                x = abs(f1(iclu,1)-f2(jclu,1));
                                z = abs(f1(iclu,2)-f2(jclu,2));
%                                 if x<=2 && z<=2
                                    C(iclu,jclu) = sqrt((f1(iclu,1)-f2(jclu,1))^2 + (f1(iclu,2)-f2(jclu,2))^2 + (f1(iclu,3)-f2(jclu,3))^2);
%                                 else
%                                     C(iclu,jclu) = 0;
%                                 end
                            end
                        end
                        C = C(:); % vectorized vertically
                        tau = 1; %partial parameter
                        b = [ones(size(f1,1),1); ones(size(f2,1),1); round(tau*min(size(f1,1), size(f2,1)))]; % parameter that A computed: all S1, all S2, smaller of total weights between S1 and S2
                        opts.nonnegativity = true;                                                 % Make sure F is non-negative
                        opts.printStopCrit = true;
                        %opts.restart =
                        %opts.maxCounts = [ Inf, Inf, iteration, Inf, Inf ];
                        ind = find(C <= inf);
                        C_new = C(ind);
                        A = @(X,T) emdConstraints_threshold(X, ind, size(f1,1),size(f2,1),T); % Constraints matrix, T is option, takes value 0, 1, and 2
                        res  = solver_sLP_Rplus(C_new, A, b, mu, [], [], opts); % Run the Generic linear programming in standard form. c and b are vectors, A is matrix

                        cost = sum(res.*C_new);                                                        % EMD value
                        P = zeros(size(f1,1),size(f2,1));
                        P(ind) = res;



                        EMD_cost(ich,i) = cost;
                        Pmatrix{ich,i} = P;
                end
                %                 % calculate norm
                %                 C = reshape(C,[length(f1),length(f2)]);
                %                 %distance(ip,i) = mean(C,'all'); %calculate mean interneuron distance
                %                 sortC = sort(C);
                %                 neighborAve(ich,i) = mean(sortC(2,:)); %find mean of all nearest neighbor distance
                %                 move = diag(C);
                %                 offsetMean(ich,i) = mean(move); %calculate mean offset between points

                noMatch = [];
                switch pars.changeType
                    case 'gain'
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


                        % find nonzero matches
                        matched_ind = [];
                        P_nz = find(sum(P,1) ~= 0); %find cols with entries
                        [~,max_ind] = max(P(:,P_nz),[],1); %max in every nonzero col
                        matched_ind(:,1) = P_nz; %f2 index has match
                        matched_ind(:,2) = max_ind; %corresponding f1 index match
                        % find cols with no match
                        noMatch(:,1) = find(sum(P,1) == 0); %f2 cols with all zeros
                        % find f1 correspondence in f2
                        f2_same = points.f2_same; %points preserved in f2
                        [q, f2_unchanged] = ismember(f2_same, f2, 'rows'); %the new index in f2 corresponds to true f1 pair index
                        f2_unchanged(length(f1)+1:length(f2),1) = 0; %added points have no pairs
                        matched_ind(:,3) = f2_unchanged(P_nz); %correct f1

                        % count matches
                        numIncorrect(ich, i) = sum(matched_ind(:,3) ~= matched_ind(:,2));
                        perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1)+ich);
                        if ~isempty(noMatch) %add number of correct unmatch
                            incorrect_no = length(noMatch) - (sum(ismember(noMatch,loss)) + sum(ismember([length(f1)+1:1:length(f1)+ich],noMatch)));
                            numIncorrect(ich, i) = numIncorrect(ich, i) + incorrect_no;
                            perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1)+ich);
                        end

                        %                         correct_dis = zeros(length(matched_ind),1);
                        %                         for ir = 1:length(matched_ind)
                        %                             if matched_ind(ir,3) == 0
                        %                                 matched_ind(ir,4) = sqrt( (f1(matched_ind(ir,2),1)-f2(matched_ind(ir,1),1))^2 + (f1(matched_ind(ir,2),2)-f2(matched_ind(ir,1),2))^2 + (f1(matched_ind(ir,2),3)-f2(matched_ind(ir,1),3))^2 );
                        %                                 correct_dis(ir,1) = 0;
                        %                             elseif matched_ind(ir,2) > length(f1)
                        %                                 matched_ind(ir,4) = 0;
                        %                                 correct_dis(ir,1) = sqrt( (f1(matched_ind(ir,3),1)-f2(matched_ind(ir,1),1))^2 + (f1(matched_ind(ir,3),2)-f2(matched_ind(ir,1),2))^2 + (f1(matched_ind(ir,3),3)-f2(matched_ind(ir,1),3))^2 );
                        %                             else
                        %                                 matched_ind(ir,4) = sqrt( (f1(matched_ind(ir,2),1)-f2(matched_ind(ir,1),1))^2 + (f1(matched_ind(ir,2),2)-f2(matched_ind(ir,1),2))^2 + (f1(matched_ind(ir,2),3)-f2(matched_ind(ir,1),3))^2 );
                        %                                 correct_dis(ir,1) = sqrt( (f1(matched_ind(ir,3),1)-f2(matched_ind(ir,1),1))^2 + (f1(matched_ind(ir,3),2)-f2(matched_ind(ir,1),2))^2 + (f1(matched_ind(ir,3),3)-f2(matched_ind(ir,1),3))^2 );
                        %                             end
                        %                         end
                        %                         figure(i)
                        %                         histogram(matched_ind(:,4),10)
                        %                         figure(10+i)
                        %                         histogram(correct_dis)



                    case 'loss'
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



                        % find nonzero matches
                        matched_ind = [];
                        noMatch = [];
                        P_nz = find(sum(P,1) ~= 0); %find cols with entries
                        [~,max_ind] = max(P(:,P_nz),[],1); %max in every nonzero col
                        matched_ind(:,1) = P_nz; %f2 index has match
                        matched_ind(:,2) = max_ind; %corresponding f1 index match
                        % find cols with no match
                        noMatch(:,1) = find(sum(P,1) == 0); %f2 cols with all zeros
                        % find f1 correspondence in f2
                        f2_same = points.f2_same; %points preserved in f2
                        [q, f2_unchanged] = ismember(f2_same, f2, 'rows'); %the new index in f2 corresponds to true f1 pair index, 0 if no pair
                        matched_ind(:,3) = f2_unchanged(P_nz); %correct f1

                        % count matches
                        numIncorrect(ich, i) = sum(matched_ind(:,3) ~= matched_ind(:,2));
                        perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1));
                        if ~isempty(noMatch) %add number of correct unmatch
                            incorrect_no = length(noMatch) - (sum(ismember(noMatch,loss)) + sum(ismember([length(f1)+1:1:length(f1)+ich],noMatch)));
                            numIncorrect(ich, i) = numIncorrect(ich, i) + incorrect_no;
                            perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1));
                        end

                        %                         correct_dis = zeros(length(matched_ind),1);
                        %                         for ir = 1:length(matched_ind)
                        %                             if matched_ind(ir,3) == 0
                        %                                 matched_ind(ir,4) = sqrt( (f1(matched_ind(ir,2),1)-f2(matched_ind(ir,1),1))^2 + (f1(matched_ind(ir,2),2)-f2(matched_ind(ir,1),2))^2 + (f1(matched_ind(ir,2),3)-f2(matched_ind(ir,1),3))^2 );
                        %                                 correct_dis(ir,1) = 0;
                        %                             elseif matched_ind(ir,2) > length(f1)
                        %                                 matched_ind(ir,4) = 0;
                        %                                 correct_dis(ir,1) = sqrt( (f1(matched_ind(ir,3),1)-f2(matched_ind(ir,1),1))^2 + (f1(matched_ind(ir,3),2)-f2(matched_ind(ir,1),2))^2 + (f1(matched_ind(ir,3),3)-f2(matched_ind(ir,1),3))^2 );
                        %                             else
                        %                                 matched_ind(ir,4) = sqrt( (f1(matched_ind(ir,2),1)-f2(matched_ind(ir,1),1))^2 + (f1(matched_ind(ir,2),2)-f2(matched_ind(ir,1),2))^2 + (f1(matched_ind(ir,2),3)-f2(matched_ind(ir,1),3))^2 );
                        %                                 correct_dis(ir,1) = sqrt( (f1(matched_ind(ir,3),1)-f2(matched_ind(ir,1),1))^2 + (f1(matched_ind(ir,3),2)-f2(matched_ind(ir,1),2))^2 + (f1(matched_ind(ir,3),3)-f2(matched_ind(ir,1),3))^2 );
                        %                             end
                        %                         end
                        %                         figure(i)
                        %                         histogram(matched_ind(:,4),10)
                        %                         figure(10+i)
                        %                         histogram(correct_dis)






                    case 'mixed'
                        % find nonzero matches
                        matched_ind = [];
                        noMatch = [];
                        P_nz = find(sum(P,1) ~= 0); %find cols with entries
                        [~,max_ind] = max(P(:,P_nz),[],1); %max in every nonzero col
                        matched_ind(:,1) = P_nz; %f2 index has match
                        matched_ind(:,2) = max_ind; %corresponding f1 index match
                        % find cols with no match
                        noMatch(:,1) = find(sum(P,1) == 0); %f2 cols with all zeros
                        % find f1 correspondence in f2
                        f2_same = points.f2_same; %points preserved in f2
                        [q, f2_unchanged] = ismember(f2_same, f2, 'rows'); %the new index in f2 corresponds to true f1 pair index, 0 if no pair
                        matched_ind(:,3) = f2_unchanged(P_nz); %correct f1

                        % count matches
                        numIncorrect(ich, i) = sum(matched_ind(:,3) ~= matched_ind(:,2));
                        perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1)+ich);
                        if ~isempty(noMatch) %add number of correct unmatch
                            incorrect_no = length(noMatch) - (sum(ismember(noMatch,loss)) + sum(ismember([length(f1)+1:1:length(f1)+ich],noMatch)));
                            numIncorrect(ich, i) = numIncorrect(ich, i) + incorrect_no;
                            perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1)+ich);
                        end

%                         correct_dis = zeros(length(matched_ind),1);
%                         for ir = 1:length(matched_ind)
%                             if matched_ind(ir,3) == 0
%                                 matched_ind(ir,4) = sqrt( (f1(matched_ind(ir,2),1)-f2(matched_ind(ir,1),1))^2 + (f1(matched_ind(ir,2),2)-f2(matched_ind(ir,1),2))^2 + (f1(matched_ind(ir,2),3)-f2(matched_ind(ir,1),3))^2 );
%                                 correct_dis(ir,1) = 0;
%                             elseif matched_ind(ir,2) > length(f1)
%                                 matched_ind(ir,4) = 0;
%                                 correct_dis(ir,1) = sqrt( (f1(matched_ind(ir,3),1)-f2(matched_ind(ir,1),1))^2 + (f1(matched_ind(ir,3),2)-f2(matched_ind(ir,1),2))^2 + (f1(matched_ind(ir,3),3)-f2(matched_ind(ir,1),3))^2 );
%                             else
%                                 matched_ind(ir,4) = sqrt( (f1(matched_ind(ir,2),1)-f2(matched_ind(ir,1),1))^2 + (f1(matched_ind(ir,2),2)-f2(matched_ind(ir,1),2))^2 + (f1(matched_ind(ir,2),3)-f2(matched_ind(ir,1),3))^2 );
%                                 correct_dis(ir,1) = sqrt( (f1(matched_ind(ir,3),1)-f2(matched_ind(ir,1),1))^2 + (f1(matched_ind(ir,3),2)-f2(matched_ind(ir,1),2))^2 + (f1(matched_ind(ir,3),3)-f2(matched_ind(ir,1),3))^2 );
%                             end
%                         end
%                         figure(i)
%                         histogram(matched_ind(:,4),10)
%                         figure(10+i)
%                         histogram(correct_dis)
                end
            end



        end
end

% numIncorrect_tfocs_mixed_revised = numIncorrect;
% save 'C:\Users\labadmin\Desktop\EMD_test_simulator\Algorithm comparison\mixed\tfocs_mixed50_revised.mat' numIncorrect_tfocs_mixed_revised


% ----------------------- Helper functions ----------------------------------------
% Generate point data
function data_generator(pars, name)
for ich = 1:length(pars.change)
    for i = 1:pars.nTrials % trials
        iseed = i;
        rng(iseed);

        f1 = zeros([pars.npts,pars.dimension]);

        f1(:,1) = pars.xRange(1) + (pars.xRange(2)-pars.xRange(1))*rand([pars.npts,1]); % random x location
        f1(:,2) = pars.zRange(1) + (pars.zRange(2)-pars.zRange(1))*rand([pars.npts,1]);
        f1(:,3) = pars.yRange(1) + (pars.yRange(2)-pars.yRange(1))*rand([pars.npts,1]);

        f2 = f1;
        f2(:,2) = f2(:,2) + pars.y_drift;
        f2(:,1) = f2(:,1) + pars.x_err_std*randn(length(f1),1); % add error: 10 horizontal + 1-5 vertical
        f2(:,2) = f2(:,2) + pars.y_err_std*randn(length(f1),1);
        f2(:,3) = f2(:,3) + pars.z_err_std*randn(length(f1),1);

        switch pars.changeType
            case 'gain'
                loss = 0;
                f2_more = [];
                f2_more(:,1) = pars.xRange(1) + (pars.xRange(2)-pars.xRange(1))*rand([pars.change(ich),1]); % random x location
                f2_more(:,2) = pars.zRange(1) + (pars.zRange(2)-pars.zRange(1))*rand([pars.change(ich),1]);
                f2_more(:,3) = pars.yRange(1) + (pars.yRange(2)-pars.yRange(1))*rand([pars.change(ich),1]);
                f2_same = f2;
                f2 = [f2; f2_more];
            case 'loss'
                f2_more = 0;
                loss = randperm(pars.npts,pars.change(ich)); %unique index of loss points
                f2_same = f2;
                f2(loss,:) = [];
                f2_same(loss,:) = NaN;
            case 'mixed'
                loss = randperm(pars.npts,pars.change(ich)); %unique of loss points
                f2_same = f2;
                f2(loss,:) = [];
                f2_same(loss,:) = NaN;
                f2_more = [];
                f2_more(:,1) = pars.xRange(1) + (pars.xRange(2)-pars.xRange(1))*rand([pars.change(ich),1]); % random x location
                f2_more(:,2) = pars.zRange(1) + (pars.zRange(2)-pars.zRange(1))*rand([pars.change(ich),1]);
                f2_more(:,3) = pars.yRange(1) + (pars.yRange(2)-pars.yRange(1))*rand([pars.change(ich),1]);
                f2 = [f2; f2_more];
        end
        %save ([name,'trial', num2str(i), ',', ' drift=', num2str(pars.y_drift), ' error=', num2str(pars.x_err_std), '&', num2str(pars.y_err_std), '&', num2str(pars.z_err_std), '.mat'], 'f1', 'f2', 'f2_same', 'f2_more', 'loss')

        idx = find(~isnan(f2_same(:,1)));
        figure(i);
        co = 1:length(f1);
        co2 = 1:length(f2);
        colormap("colorcube")
        scatter3(f1(:,1),f1(:,2),f1(:,3),50,co,'o'); hold on; % plot original set
        scatter3(f2(:,1),f2(:,2),f2(:,3),20,co2,'filled','d'); hold on; % plot drifted set
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        flowfig = quiver3(f1(idx,1), f1(idx,2), f1(idx,3), f2_same(idx,1)-f1(idx,1), f2_same(idx,2)-f1(idx,2), f2_same(idx,3)-f1(idx,3), 'off'); hold off;
        %savefig([name,'trial', num2str(i), ',', ' drift=', num2str(pars.y_drift), ' error=', num2str(pars.x_err_std), '&', num2str(pars.y_err_std), '&', num2str(pars.z_err_std), '.fig'])
    end
end
end