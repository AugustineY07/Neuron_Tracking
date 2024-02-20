function [centCorr, centL2] = calcCenteredCorrL2( mw1, mw2, chan_pos, zStep, nRow )

% simple version -- find center row for each unit, get -nRow:centerRow:nRow
% sites to include in correlation.

% Note that we need the channel positions, it is assumed that mw1 and mw2
% have identical positions

% recommend nRow = 5
% for edge cases, build the same size neighborhood for the two units

pp_all_1 = squeeze(max(mw1,[],2)-min(mw1,[],2)); %for 2 units, find peak site
pp_all_2 = squeeze(max(mw2,[],2)-min(mw2,[],2));
[~, pk_chan_1] = max(pp_all_1); %ptp, peak site
[~, pk_chan_2] = max(pp_all_2);

% probe geometry
minZ = min(chan_pos(:,2)); %min z of current probe geometry
maxZ = max(chan_pos(:,2)); %max z
%minX_pos = min(chan_pos(:,1));

unitZ = [chan_pos(pk_chan_1,2), chan_pos(pk_chan_2,2)]; %z location of the two peak channels
%dZ = unitZ(2)-unitZ(1); %difference

% for each unit, calculate the channel z values that will be included
zIncl = calcRows(unitZ,minZ,maxZ,zStep,nRow);
% fprintf('zIncl = %d\n', zIncl);


cv = 0;         %for summing correlations
cl2 = 0;        %for summing L2 norm
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
            wf1 = mw1(cc1(j),:);
            wf2 = mw2(cc2(j),:);
            denom = max(norm(wf1), norm(wf2));
            if denom > 0
                nSite = nSite + 1;
                cl2 = cl2 + norm(wf2-wf1)/max(norm(wf1), norm(wf2));
                % commenting out correlation; not currently used
                % cv = cv + corr(wf1',wf2');  
            end
        end
    end
end
% fprintf('cv = %d\n', cv);
centCorr = cv/nSite;
% fprintf('nSite = %d\n', nSite);
centL2 = cl2/nSite;
% fprintf('centL2 = %d\n', centL2);
end



function zIncl = calcRows(unitZ,minZ,maxZ,zStep,neighRow)

if min(unitZ)-neighRow*zStep >= minZ && max(unitZ)+neighRow*zStep <= maxZ
    % calc region is in the interior of probe
    nRow = 2*neighRow+1;
    zIncl = zeros(2,nRow);
    zIncl(1,:) = unitZ(1)-neighRow*zStep + (0:nRow-1)*zStep;
    zIncl(2,:) = unitZ(2)-neighRow*zStep + (0:nRow-1)*zStep;

elseif (min(unitZ)-neighRow*zStep < minZ) && (max(unitZ)+neighRow*zStep <= maxZ)
    % Region is clipped at the bottom for at least one unit.
    % Neither is clipped at the top.
    % the two will be set to include the same number of rows
    rowsBelow = (min(unitZ)-minZ)/zStep;  % number of rows available below the peak for the clip
    nRow = neighRow+1+rowsBelow;
    zIncl = zeros(2,nRow);
    zIncl(1,:) = unitZ(1)-rowsBelow*zStep + (0:nRow-1)*zStep;
    zIncl(2,:) = unitZ(2)-rowsBelow*zStep + (0:nRow-1)*zStep;

elseif (max(unitZ)+neighRow*zStep > maxZ) && (min(unitZ)-neighRow*zStep >= minZ)
    % Region is clipped at the top for at least one unit.
    % Neither is clipped at the bottom.
    % the two will be set to include the same number of rows
    rowsAbove = (maxZ-max(unitZ))/zStep;
    nRow = neighRow+1+rowsAbove;
    zIncl = zeros(2,nRow);
    zIncl(1,:) = unitZ(1)-neighRow*zStep + (0:nRow-1)*zStep;
    zIncl(2,:) = unitZ(2)-neighRow*zStep + (0:nRow-1)*zStep;

elseif (min(unitZ)-neighRow*zStep < minZ) && (max(unitZ)+neighRow*zStep > maxZ)
    % Region is clipped at the top for one unit, at the bottom for the
    % other unit
    rowsAbove = (maxZ-max(unitZ))/zStep;
    rowsBelow = (min(unitZ)-minZ)/zStep;
    nRow = 1+rowsBelow+rowsAbove;
    zIncl(1,:) = unitZ(1)-rowsBelow*zStep + (0:nRow-1)*zStep;
    zIncl(2,:) = unitZ(2)-rowsBelow*zStep + (0:nRow-1)*zStep;
end


end



