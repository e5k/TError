% TError package
% Power-law fit for volume calulation with the method of Bonadonna and Houghton 2005
% Note that we use here the commonly used approximation of an exponential
% fit of log10(x), log10(y) to approximate a power-law
% xdata: Square root of the area (km)
% ydata: Thickness (cm)
% T0:    Intercept obtained with the method of Fierstein and Nathenson (1992)
% C:     Distal integration limit (km)
function [vol, fit_BH05_v] = bh2005(xdata, ydata, T0, C)

N   = size(xdata,1);

xdata = log10(xdata);
ydata = log10(ydata);

fit_BH05_v(2) = ( sum(xdata.*ydata) - sum(xdata)*sum(ydata)/N ) / ( sum(xdata.^2) - sum(xdata)^2/N ); % Slope
fit_BH05_v(1) = mean(ydata)-fit_BH05_v(2)*mean(xdata);

vol                 = get_vol_BH05(exp(T0), -1*fit_BH05_v(2), 10^fit_BH05_v(1), C);