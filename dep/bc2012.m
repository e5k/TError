
% TError package
% Weibull fit for volume calulation with the method of Bonadonna and Costa 2012
% xdata: Square root of the area (km)
% ydata: Thickness (cm)
% lam_r: Range of lambda values used in optimisation algorithm
% n_r:   Range of n values used in optimisation algorithm
function [vol, fit_BC12_v] = bc2012(xdata, ydata, lam_r, n_r)

wbl                 = @(x02)weibe(x02, xdata, ydata); % Weibull function
[fit_BC12_v]        = fminsearchbnd(wbl, [.5, .5, .5], [.1, lam_r(1) ,n_r(1)], [5000, lam_r(2), n_r(2)]); % Optimization algorithm

fit_BC12_v(1)       = fit_BC12_v(2)^(fit_BC12_v(3)-2) * sum( (xdata.^(fit_BC12_v(3)-2) ./ ydata).*exp(-(xdata/fit_BC12_v(2)).^fit_BC12_v(3)) ) .* ...
    (sum( ((xdata.^(fit_BC12_v(3)-2) ./ ydata).*exp(-(xdata./fit_BC12_v(2)).^fit_BC12_v(3)) ).^2)).^-1; % Theta calculated from eq 2 of Daggit et al.

vol                 = get_vol_BC12(fit_BC12_v(1), fit_BC12_v(2), fit_BC12_v(3));

% Weilbull fit
function F = weibe(x, xdata, ydata)
x(1) = x(2)^(x(3)-2) * sum( (xdata.^(x(3)-2) ./ ydata).*exp(-(xdata/x(2)).^x(3)) ) .* ...
    (sum( ((xdata.^(x(3)-2) ./ ydata).*exp(-(xdata./x(2)).^x(3)) ).^2)).^-1; % Theta calculated from eq 2 of Daggit et al.


y  = x(1)*(xdata/x(2)).^(x(3)-2).*exp(-1*(xdata/x(2)).^x(3));
F  = sum(((ydata - y)./ydata).^2) + log(sum(((ydata - y)./ydata).^2));