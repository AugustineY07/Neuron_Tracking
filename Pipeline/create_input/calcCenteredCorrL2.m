function [centCorr, centL2] = calcCenteredCorrL2( mw1, mw2, chanPos, nRow )

% simple version -- find center row for each unit, get -nRow:centerRow:nRow
% sites to include in correlation.

% Note that we need the channel positions, it is assumed that mw1 and mw2
% have identical positions

% recommend nRow = 5
% for edge cases, build the same size neighborhood for the two units

pp_all_1 = squeeze(max(mw1,[],2)-min(mw1,[],2));
pp_all_2 = squeeze(max(mw2,[],2)-min(mw2,[],2));
[pp_unit_1, pk_chan_1] = max(pp_all_1);
[pp_unit_2, pk_chan_2] = max(pp_all_2);

% probe geometry, np 2.0 specific
minZ = min(chanPos(:,2));
maxZ = max(chanPos(:,2));
zStep = 15;
minX_pos = min(chanPos(:,1));
xStep = 32;

unitZ = [chanPos(pk_chan_1,2), chanPos(pk_chan_2,2)];
dZ = unitZ(2)-unitZ(1);

% for each unit, calculate the z values that will be included 
zIncl = calcRows(unitZ,minZ,maxZ,nRow);

corrVals = zeros([11,1]); % to accumulate corrlation values across tshifts

    cv = 0;         %for summing correlations
    cl2 = 0;        %fo summing L2 norm
    nSite = 0;
    for i = 1:numel(zIncl(1,:))
        % specific to 2 column arrangement of sites
        currZ = [zIncl(1,i), zIncl(1,i) + dZ];
        % up to two sites, for the two columns
        cc1 = find(ismember(chanPos(:,2),currZ(1)));
        cc2 = find(ismember(chanPos(:,2),currZ(2)));  
        cc2_xPos = chanPos(cc2,:);
        for j = 1:numel(cc1)
            cX = chanPos(cc1(j));
            if ismember(cX, cc2_xPos)
                nSite = nSite + 1;
                %fprintf('chans for correlation: %d, %d\n', cc1(j), cc2(find(cc2_xPos==cX)));
                w1 = mw1(cc1(j),:);
                w2 = mw2(cc2(find(cc2_xPos==cX)),:);  
                cv = cv + corr(w1',w2');
                %fprintf('norm w1-w2, norm w1, norm w2: %.2f, %.2f, %.2f\n', norm(w1-w2), norm(w1), norm(w2));
                %cl2 = cl2 + (1 - norm(w2-w1)/max(norm(w1), norm(w2)));
                cl2 = cl2 + norm(w2-w1)/max(norm(w1), norm(w2));
            end
        end
    end
centCorr = cv/nSite;
centL2 = cl2/nSite;
end



function zIncl = calcRows(unitZ,minZ,maxZ,neighRow)

zStep = 15;
total_rows = 1 + (maxZ - minZ)/zStep;

rowInd = (unitZ - minZ)/zStep; % zero based

if min(rowInd) >= neighRow && max(rowInd) <= total_rows - neighRow
    % calc region is in the interior of probe
    zIncl(1,:) = minZ + [rowInd(1)-neighRow:1:rowInd(1)+neighRow]*zStep;
    zIncl(2,:) = minZ + [rowInd(2)-neighRow:1:rowInd(2)+neighRow]*zStep;
elseif min(rowInd) < neighRow
    % region is clipped at the bottom for at least one unit
    % the two will be set to include the same number of rows
    rowsBelow = min(rowInd);
    zIncl(1,:) = minZ + [rowInd(1)-rowsBelow:1:rowInd(1)+neighRow]*zStep;
    zIncl(2,:) = minZ + [rowInd(2)-rowsBelow:1:rowInd(2)+neighRow]*zStep;
elseif max(rowInd) > total_rows-neighRow
    % region is clipped at the top for at least one unit
    % the two will be set to include the same number of rows
    rowsAbove = total_rows - max(rowInd);
    zIncl(1,:) = minZ + [rowInd(1)-neighRow:rowInd(1)+rowsAbove]*zStep;
    zIncl(2,:) = minZ + [rowInd(2)-neighRow:rowInd(2)+rowsAbove]*zStep;
end


end



