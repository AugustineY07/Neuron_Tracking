function [centCorr, centL2] = calcCenteredCorrL2( mw1, mw2, chan_pos, zStep, nRow )

% simple version -- find center row for each unit, get -nRow:centerRow:nRow
% sites to include in correlation.

% Note that we need the channel positions, it is assumed that mw1 and mw2
% have identical positions

% recommend nRow = 5
% for edge cases, build the same size neighborhood for the two units

pp_all_1 = squeeze(max(mw1,[],2)-min(mw1,[],2)); %for 2 units, find peak site
pp_all_2 = squeeze(max(mw2,[],2)-min(mw2,[],2));
[pp_unit_1, pk_chan_1] = max(pp_all_1); %ptp, peak site
[pp_unit_2, pk_chan_2] = max(pp_all_2);

% probe geometry, (np 2.0 specific,) geometry non-specific
minZ = min(chan_pos(:,2)); %min z of current probe geometry
maxZ = max(chan_pos(:,2)); %max z
minX_pos = min(chan_pos(:,1));

unitZ = [chan_pos(pk_chan_1,2), chan_pos(pk_chan_2,2)]; %z location of the two peak channels
dZ = unitZ(2)-unitZ(1); %difference

% for each unit, calculate the channel z values that will be included
zIncl = calcRows(chan_pos,unitZ,minZ,maxZ,zStep,nRow);
fprintf('zIncl = %d\n', zIncl);


cv = 0;         %for summing correlations
cl2 = 0;        %fo summing L2 norm
nSite = 0;
mid_point = (max(chan_pos(:,1))-min(chan_pos(:,1))) / 2 + min(chan_pos(:,1));
for i = 1:numel(zIncl(1,:)) %row loop = 11
    % (specific to 2 column arrangement of sites)
    %     currZ = [zIncl(1,i), zIncl(1,i) + dZ]; %z of each channel included
    currZ = [zIncl(1,i), zIncl(2,i)]; %z of each channel included
    % (up to two sites, for the two columns)
    cc1 = find(ismember(chan_pos(:,2),currZ(1))); %channel index among all chan_pos of unit 1
    cc2 = find(ismember(chan_pos(:,2),currZ(2))); %channel index among all chan_pos of unit 2
    %         cc2_xPos = chan_pos(cc2,:); %the two sites with same z of unit 2
    for j = 1:min(numel(cc1),numel(cc2)) %col loop = 2, some rows may have only 1 site
        cX1 = chan_pos(cc1(j));
        cX2 = chan_pos(cc2(j));
        %corr of matching sites
        if (cX1 < mid_point && cX2 < mid_point) || (cX1 > mid_point && cX2 > mid_point)%ismember(cX, cc2_xPos)
            nSite = nSite + 1;
            %fprintf('chans for correlation: %d, %d\n', cc1(j), cc2(find(cc2_xPos==cX)));
            wf1 = mw1(cc1(j),:);
            wf2 = mw2(cc2(j),:);
            %                 wf2 = mw2(cc2(find(cc2_xPos==cX)),:);
            cv = cv + corr(wf1',wf2');
            %fprintf('norm wf1-wf2, norm wf1, norm wf2: %.2f, %.2f, %.2f\n', norm(wf1-wf2), norm(wf1), norm(wf2));
            %cl2 = cl2 + (1 - norm(wf2-wf1)/max(norm(wf1), norm(wf2)));
            cl2 = cl2 + norm(wf2-wf1)/max(norm(wf1), norm(wf2));
        end
    end
end
fprintf('cv = %d\n', cv);
centCorr = cv/nSite;
fprintf('nSite = %d\n', nSite);
centL2 = cl2/nSite;
fprintf('centL2 = %d\n', centL2);
end




%neighRow=5
function zIncl = calcRows(chan_pos,unitZ,minZ,maxZ,zStep,neighRow)

total_rows = 1 + (maxZ - minZ)/zStep; %number of rows in the probe

% rowInd = (unitZ - minZ)/zStep; % zero based
% find the row of pk channels when chan_map not always in order
rowInd = [find(chan_pos(:,2) == unitZ(1),1) find(chan_pos(:,2) == unitZ(2),1)];


if min(unitZ)-neighRow*zStep >= minZ && max(unitZ)+neighRow*zStep <= maxZ
    for ir = 1:2*neighRow+1
        % calc region is in the interior of probe
        zIncl(1,ir) = unitZ(1)-neighRow*zStep + (ir-1)*zStep;
        zIncl(2,ir) = unitZ(2)-neighRow*zStep + (ir-1)*zStep;
    end
elseif (min(unitZ)-neighRow*zStep < minZ) && (max(unitZ)+neighRow*zStep <= maxZ)
    % one region is clipped at the top for at least one unit, but the other
    % is not clipped at the bottom
    % the two will be set to include the same number of rows
    rowsAbove = (min(unitZ)-minZ)/zStep;
    for ir = 1:neighRow+1+rowsAbove
        % calc region is in the interior of probe
        zIncl(1,ir) = unitZ(1)-rowsAbove*zStep + (ir-1)*zStep;
        zIncl(2,ir) = unitZ(2)-rowsAbove*zStep + (ir-1)*zStep;
    end
elseif (max(unitZ)+neighRow*zStep > maxZ) && (min(unitZ)-neighRow*zStep >= minZ)
    % one region is clipped at the bottom for at least one unit, the other
    % is not clipped at the top
    % the two will be set to include the same number of rows
    rowsBelow = (maxZ-max(unitZ))/zStep;
    for ir = 1:neighRow+1+rowsBelow
        % calc region is in the interior of probe
        zIncl(1,ir) = unitZ(1)-neighRow*zStep + (ir-1)*zStep;
        zIncl(2,ir) = unitZ(2)-neighRow*zStep + (ir-1)*zStep;
    end
elseif (min(unitZ)-neighRow*zStep < minZ) && (max(unitZ)+neighRow*zStep > maxZ)
    % one region cipped at top and other at bottom
    rowsBelow = (maxZ-max(unitZ))/zStep;
    rowsAbove = (min(unitZ)-minZ)/zStep;
    for ir = 1:1+rowsBelow+rowsAbove
        zIncl(1,ir) = unitZ(1)-rowsAbove*zStep + (ir-1)*zStep;
        zIncl(2,ir) = unitZ(2)-rowsAbove*zStep + (ir-1)*zStep;
    end
end



% if min(rowInd) >= neighRow && max(rowInd) <= total_rows - neighRow
%     % calc region is in the interior of probe
%     zIncl(1,:) = minZ + [rowInd(1)-neighRow:1:rowInd(1)+neighRow]*zStep;
%     zIncl(2,:) = minZ + [rowInd(2)-neighRow:1:rowInd(2)+neighRow]*zStep;
% elseif min(rowInd) < neighRow
%     % region is clipped at the bottom for at least one unit
%     % the two will be set to include the same number of rows
%     rowsBelow = min(rowInd);
%     zIncl(1,:) = minZ + [rowInd(1)-rowsBelow:1:rowInd(1)+neighRow]*zStep;
%     zIncl(2,:) = minZ + [rowInd(2)-rowsBelow:1:rowInd(2)+neighRow]*zStep;
% elseif max(rowInd) > total_rows-neighRow
%     % region is clipped at the top for at least one unit
%     % the two will be set to include the same number of rows
%     rowsAbove = total_rows - max(rowInd);
%     zIncl(1,:) = minZ + [rowInd(1)-neighRow:rowInd(1)+rowsAbove]*zStep;
%     zIncl(2,:) = minZ + [rowInd(2)-neighRow:rowInd(2)+rowsAbove]*zStep;
% end


end



