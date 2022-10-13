% This function contains all EMD calculation with all options

clear all;
close all;

% CHANGE this to your path
addpath(genpath('C:\Users\labadmin\Desktop\EMD_test_simulator')) % path to simulator folder
addpath(genpath('D:\Data\Pipeline\EMD\Localization_EMD_TFOCS')) 
addpath(genpath('D:\Data\Pipeline\EMD\TFOCS\TFOCS-master')) % TFOCS package

pars.nTrials = 10;
pars.npts = 50; %number of points per trial
pars.xRange = [-50,50];
pars.zRange = [-50,770];
pars.yRange = [-10,10];
pars.y_drift = 0;
pars.pos_err_std = 5;
pars.x_err_std = 1;
pars.y_err_std = 1;
pars.z_err_std = 1;
pars.plot = 1;
pars.save = 1; %0 = not save, 1 = save data

pars.rootD = 'C:\Users\labadmin\Desktop\EMD_test_simulator\'; % the simulated data files are in this folder, CHANGE to your path
pars.runType = 1; %0 = data generation, 1 = run EMD
pars.alg = 'TFOCS'; %1 = matlab, 2 = CVX, 3 = TFOCS
pars.dimension = 3; %2 = '2D', 3 = '3D'
pars.errorType = 'xyz'; %'Pos' = pos_err_std, 'xyz' = x_y_z_error
pars.mode = 'gainloss'; %0 = standard, 1 = lumpiness, 2 = gainloss
pars.block = 4; %lumpiness: number of blocks
pars.changeType = 'loss'; %gainloss: 1 = gain, 2 = loss, 3 = mixed
pars.change = [0.1]; %number of lost/gain points
pars.lost = 0.1;

% CHANGE to your path name
switch pars.mode
    case {'standard','lumpiness'}
        name = ['C:\Users\labadmin\Desktop\EMD_test_simulator\data\lost0.2gain0.1\', num2str(pars.dimension),'\', pars.errorType,'\', pars.mode,'\'];
    case 'gainloss'
        name = ['C:\Users\labadmin\Desktop\EMD_test_simulator\data\', num2str(pars.dimension),'\', pars.errorType,'\', pars.mode,'\', pars.changeType,'\'];
        %lost0.2gain0.1\
end

switch pars.runType %data v.s run
    case 0
        data_generator(pars, name);
    case 1
        numIncorrect = [];
        perIncorrect = [];

        for ich = 1%1:length(pars.change)
            for i = 1:pars.nTrials
                % 100 points dataset
                %points = load([name,'trial', num2str(i), ',', ' drift=', num2str(pars.y_drift), ' error=', num2str(pars.x_err_std), '&', num2str(pars.y_err_std), '&', num2str(pars.z_err_std), '.mat']);
                % 50 points dataset
                points = load([pars.rootD, pars.mode, '\data\', pars.changeType, '\trial', num2str(i), ',', 'npts', num2str(pars.npts), pars.changeType, num2str(5), ',error=', num2str(pars.x_err_std), '&', num2str(pars.y_err_std),'&', num2str(pars.z_err_std), '.mat']);

                f1 = points.f1;
                f2 = points.f2;
                f2_more = points.f2_more; %num added points
                loss = points.loss; %f1 idx of lost points
                f2_same = points.f2_same; %f2 with changed points marked
                idx_same = find(~isnan(f2_same(:,1))); %idx of unchanged points
                if loss ~= 0
                    f1_same = f1;
                    f1_same(loss,:) = NaN;
                end

                % set weights = 1
                w1 = ones(length(f1),1);%ones([length(f1),1])/length(f1);
                w2 = ones(length(f2),1);%ones([length(f2),1])/length(f2);

                % initialize distance mat
                C = [];

                switch pars.alg
                    case 'matlab'
                        [x, fval] = emd(f1, f2, w1, w2, @gdf);
                        P = reshape(x,[length(f2),length(f1)]);
                        C = gdm(f1, f2, @gdf); % Distance matrix

                        cost = sum(C.*x);
                        EMD_cost(ich,i) = cost;
                        Pmatrix{ich,i} = P;
                    case 'CVX'
                        C = gdm(f1, f2, @gdf);

                        % Solve
                        tic
                        cvx_begin quiet
                        % cvx_solver
                        cvx_solver SDPT3 %optimization solver, if slow, try smaller x first
                        cvx_precision high
                        variable P(length(f2),length(f1)) nonnegative; %f_ij >= 0
                        minimize(C(:)'*P(:)); %vectorize C by col, C became f2*f1
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

                    case 'TFOCS'

                        [P,cost] = runTfocsMatching(f1,f2,'n');
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

                matched_ind = [];
                noMatch = [];
                switch pars.changeType
                    case 'gain'
                        % find nonzero matches
                        P_nz = find(sum(P,1) ~= 0); %find rows with entries
                        [~,max_ind] = max(P(:,P_nz),[],1); %max in every nonzero row
                        matched_ind(:,1) = P_nz; %f2 index has match
                        matched_ind(:,2) = max_ind; %corresponding f1 index match
                        % find cols with no match
                        noMatch(:,1) = find(sum(P,1) == 0); %f2 rows with all zeros
                        % find f1 correspondence in f2
                        [q,  f2_unchanged] = ismember(f2_same, f2, 'rows'); %the new index in f2 corresponds to true f1 pair index
                        f2_unchanged(length(f1)+1:length(f2),1) = 0; %added points have no pairs
                        matched_ind(:,3) = f2_unchanged(P_nz); %correct f1

                        % count matches
                        incorrect_no = length(noMatch) - sum(ismember([length(f1)+1:1:length(f1)+ich],noMatch));
                        numIncorrect(ich, i) = sum(matched_ind(:,2) ~= matched_ind(:,3));
                        perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1)+ich);



                    case 'loss'
                        % find nonzero matches
                        P_nz = find(sum(P,1) ~= 0); %find f1 index with entries
                        [~,max_ind] = max(P(:,P_nz),[],1); %max in every nonzero f1 index
                        matched_ind(:,1) = P_nz; %f1 index has match
                        matched_ind(:,2) = max_ind; %corresponding f2 index match
                        % find cols with no match
                        noMatch(:,1) = find(sum(P,1) == 0); %f1 index with no match
                        % find f1 correspondence in f2
                        [q, f2_unchanged] = ismember(f2_same, f2, 'rows'); %the new index in f2 corresponds to true f1 pair index, 0 if no match
                        matched_ind(:,3) = f2_unchanged(P_nz); %correct f2 match of the nonzero f1 index

                        % count matches
                        numIncorrect(ich, i) = sum(matched_ind(:,2) ~= matched_ind(:,3));
                        perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1));
                        if ~isempty(noMatch) %add number of correct unmatch
                            incorrect_no = length(noMatch) - sum(ismember(noMatch,loss));
                            numIncorrect(ich, i) = numIncorrect(ich, i) + incorrect_no;
                            perIncorrect(ich, i) = numIncorrect(ich, i)/(length(f1));
                        end


                    case 'mixed'
                        % find nonzero matches
                        P_nz = find(sum(P,1) ~= 0); %find cols with entries
                        [~,max_ind] = max(P(:,P_nz),[],1); %max in every nonzero col
                        matched_ind(:,1) = P_nz; %f2 index has match
                        matched_ind(:,2) = max_ind; %corresponding f1 index match
                        % find cols with no match
                        noMatch(:,1) = find(sum(P,1) == 0); %f2 cols with all zeros
                        % find f1 correspondence in f2
                        [q, f2_unchanged] = ismember(f2_same, f2, 'rows'); %the new index in f2 corresponds to true f1 pair index, 0 if no pair
                        if pars.change ~= pars.lost
                            f2_unchanged(length(f1)+1:length(f2),1) = 0; %added points have no pairs
                        end
                        matched_ind(:,3) = f2_unchanged(P_nz); %correct f1

                        % count matches
                        numIncorrect(ich, i) = sum(matched_ind(:,2) ~= matched_ind(:,3));
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

                % calculate gt
                [emd_match, emd_noMatch] = gt_best_dist(pars, f1, f2, loss, f2_more, f2_same, idx_same);
                emd_gt(ich, i) = (emd_match + emd_noMatch);
            end
        end
        plot_acc_gt(perIncorrect, emd_gt, EMD_cost);
        acc_mean = mean(perIncorrect)
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
                f2_more(:,1) = pars.xRange(1) + (pars.xRange(2)-pars.xRange(1))*rand([pars.change(ich)*pars.npts,1]); % random x location
                f2_more(:,2) = pars.zRange(1) + (pars.zRange(2)-pars.zRange(1))*rand([pars.change(ich)*pars.npts,1]);
                f2_more(:,3) = pars.yRange(1) + (pars.yRange(2)-pars.yRange(1))*rand([pars.change(ich)*pars.npts,1]);
                f2_same = f2;
                f2 = [f2; f2_more];
            case 'loss'
                f2_more = 0;
                loss = randperm(pars.npts,pars.change(ich)*pars.npts); %unique index of loss points
                f2_same = f2;
                f2(loss,:) = [];
                f2_same(loss,:) = NaN;
            case 'mixed'
                loss = randperm(pars.npts,pars.lost*pars.npts); %unique of loss points
                f2_same = f2;
                f2(loss,:) = [];
                f2_same(loss,:) = NaN;
                f2_more = [];
                f2_more(:,1) = pars.xRange(1) + (pars.xRange(2)-pars.xRange(1))*rand([pars.change(ich)*pars.npts,1]); % random x location
                f2_more(:,2) = pars.zRange(1) + (pars.zRange(2)-pars.zRange(1))*rand([pars.change(ich)*pars.npts,1]);
                f2_more(:,3) = pars.yRange(1) + (pars.yRange(2)-pars.yRange(1))*rand([pars.change(ich)*pars.npts,1]);
                f2 = [f2; f2_more];
        end

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
        if pars.save == 1
            if ~exist([name], 'dir')
                mkdir([name])
            end
            save ([name,'trial', num2str(i), ',', ' drift=', num2str(pars.y_drift), ' error=', num2str(pars.x_err_std), '&', num2str(pars.y_err_std), '&', num2str(pars.z_err_std), '.mat'], 'f1', 'f2', 'f2_same', 'f2_more', 'loss')
            savefig([name,'trial', num2str(i), ',', ' drift=', num2str(pars.y_drift), ' error=', num2str(pars.x_err_std), '&', num2str(pars.y_err_std), '&', num2str(pars.z_err_std), '.fig'])
        end
    end
end
end


% plot
function plot_acc_gt(perIncorrect, emd_gt, EMD_cost)
% compare gt and alg
figure;
diff(:,:) = emd_gt - EMD_cost;
scatter(perIncorrect, diff, 'filled')
end


% Compare true distance v.s. algorithm result
function [emd_match, emd_noMatch] = gt_best_dist(pars, f1, f2, loss, f2_more, f2_same, idx_same)
% calculate the sum of distances among the matched pairs = gt_emd
emd_match = 0;
emd_noMatch = 0;
f1_nNoMatch = length(loss);
f2_nNoMatch = size(f2_more,1);
for im = 1:length(idx_same)
    emd_match = emd_match + gdf(f2_same(idx_same(im),:), f1(idx_same(im),:));
end

switch pars.changeType
    case 'mixed'
        if f1_nNoMatch > 0 && f2_nNoMatch > 0
            f1_unmatched = f1(loss,:);
            f2_unmatched = f2_more;
            % get the matlab emd result from the unmatched points
            w3 = ones([f1_nNoMatch,1]);
            w4 = ones([f2_nNoMatch,1]);

            C2 = gdm(f1_unmatched, f2_unmatched, @gdf);
            matched_ind2 = [];
            noMatch2 = [];

            switch pars.alg
                case 'matlab'
                    [x2, fval2] = emd(f1_unmatched, f2_unmatched, w3, w4, @gdf);
                    P2 = reshape(x2,[length(f2_unmatched),length(f1_unmatched)]);
                case 'CVX'
                    % Solve
                    tic
                    cvx_begin quiet
                    % cvx_solver
                    cvx_solver SDPT3 %optimization solver, if slow, try smaller x first
                    cvx_precision high
                    variable P2(length(f1_unmatched),length(f2_unmatched)) nonnegative; %f_ij >= 0
                    minimize(C2(:)'*P2(:)); %vectorize C by col, C became f2*f1
                    subject to
                    sum(P2,2)  <= w3(:); %can't move more to destination
                    sum(P2,1)' <= w4(:); %can't move more from origin
                    sum(P2(:)) == min(sum(w4),sum(w3)); %total amount is the smaller mass
                    cvx_end

                    stats2.emd     =   sum(C2(:)'*P2(:));
                    stats2.runtime = toc;
                    stats2.nbr_iter = cvx_slvitr; % the number of iterations taken by the solver
                case 'TFOCS'
                    mu = 1;
                    A2 = @(X,T) emdConstraints(X,size(f1_unmatched,1),size(f2_unmatched,1),T); % Constraints matrix, T is option, takes value 0, 1, and 2
                    b2 = [ones(size(f1_unmatched,1),1); ones(size(f2_unmatched,1),1); min(size(f1_unmatched,1), size(f2_unmatched,1))]; % parameter that A computed: all S1, all S2, smaller of total weights between S1 and S2

                    opts.nonnegativity = true;                                                 % Make sure F is non-negative
                    res2 = solver_sLP(C2, A2, b2, mu, [], [], opts);                               % Run the Generic linear programming in standard form. c and b are vectors, A is matrix
                    cost2 = sum(res2.*C2);                                                        % EMD value
                    P2 = reshape(res2,size(f1_unmatched,1),size(f2_unmatched,1));                                    % Reshape output
            end
            % for each row, find the max column
            P_nz2 = find(sum(P2,1) ~= 0); %find cols with entries
            [~,max_ind2] = max(P2(:,P_nz2),[],1); %max in every nonzero col
            matched_ind2(:,1) = P_nz2; %f2 index has match
            matched_ind2(:,2) = max_ind2; %corresponding f1 index match
            % find cols with no match
            noMatch2(:,1) = find(sum(P2,1) == 0); %f2 cols with all zeros
            for i2 = 1:min([f1_nNoMatch,f2_nNoMatch,size(P_nz2,2)])
                emd_noMatch = emd_noMatch + gdf(f1_unmatched(max_ind2(i2),:), f2_unmatched(i2,:));
            end
        end
end
end

