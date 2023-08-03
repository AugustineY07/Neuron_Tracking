function fitvals = fit_loc(uI, pp_all, chan_pos)

% pp_all is the set of all peak-to-peak voltages for all sites, all units
% chan_pos is the channel positions for all sites
% Might be able to improve efficiency by moving creation of fittype outside
% the function

bPlot = 0;  % make plots of each fit
npts = 10;  % number of sites to include in the fit
vpp_unit = squeeze(pp_all(uI,:,:))';
[~, pk_chan] = max(vpp_unit);

pk_x = chan_pos(pk_chan,1);
pk_z = chan_pos(pk_chan,2);

distsq = (chan_pos(:,1)-pk_x).^2 + (chan_pos(:,2)-pk_z).^2;
[~,dist_order] = sort(distsq);

chan_sel = dist_order(1:npts); % nearest npts sites
vpp_sel = vpp_unit(chan_sel);
x_sel = chan_pos(chan_sel,1);
z_sel = chan_pos(chan_sel,2);


a0 = pk_x;    %xsoma
b0 = pk_z;  %zsoma
c0 = 20;    %ysoma
d0 = 3000;  %alpha

up_bound = [a0+100, b0+100, c0+50, inf];
low_bound = [a0-100, b0-100, 0, 0];

fitL = fittype('vpp_loc(xsite, zsite, a, b, c, d)', ...
    'independent', {'xsite', 'zsite'}, 'dependent', 'vpp' );

newfit = fit( [x_sel, z_sel], vpp_sel, fitL, 'StartPoint', [a0, b0, c0, d0], ...
    'upper', up_bound, 'lower', low_bound);

if bPlot
    vpp_calc = feval(newfit,[x_sel, z_sel]);
    figure(uI)
    scatter( x_sel, z_sel, 200, vpp_sel, 'filled', 'square');
    hold on;
    scatter( x_sel+5, z_sel, 200, vpp_calc, 'filled', 'square');
    colorbar
    hold off;
    
    figure(1000 + uI)
    residual = vpp_calc - vpp_sel;
    scatter( x_sel, z_sel, 200, residual, 'filled', 'square');
    colorbar
    caxis([-50,50]);
end

fitvals = coeffvalues(newfit);

% fprintf( 'fit coefficients: \n');
% fprintf( 'x soma: %.1f\n', fitvals(1));
% fprintf( 'z soma: %.1f\n', fitvals(2));
% fprintf( 'y soma: %.1f\n', fitvals(3));
% fprintf( 'alpha: %.1f\n', fitvals(4));

 
end