function vpp = vpp_loc(xsite, zsite, a,b,c,d)

% Julien Boussard - style fit of peak-to-peak voltage vs position
% From his python code:
%   def minimize_ls(self, vec, wfs_0, z_initial, channels):
%   return wfs_0.ptp(1)-vec[3]/(((self.geom_array[channels] - [vec[0], z_initial+vec[1]])**2).sum(1) + vec[2]**2)**0.5 # vec[0]

% In matlab, can express as a surface fit (rather than minimization)
% 'Fit objects' in MATLAB order coeffcients alphabetically, so they are
% nicknamed (a,b,c,d) for the input
% 1/R model for calculating Vpp at (x,z) of a site given x,z,y and alpha of
% a given neuron.
% params (1,2,3) = soma (X, Z, Y)
% params (4) = alpha
% f = 
% alpha/( (x_site-x_soma)^2 + (z_site - z_soma)^2 + (y_soma)^2)^(1/2)

xsoma = a;
zsoma = b;
ysoma = c;
alpha = d;

vpp = zeros(size(xsite));
for i = 1:numel(xsite)
    vpp(i) = alpha/(((xsite(i) - xsoma).^2 + (zsite(i) - zsoma).^2 + ysoma.^2).^0.5);
end

end