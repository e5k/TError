% TError package
% This function plots all individual fits performed during Monte Carlo
% simulations
% xdata:    Square root of the area (km)
% Aip:      Break-in-slopes as indices
% C:        Distal integration limit (km)
% fit_FN92: Exponential fit
% fit_BH05: Power-law fit
% fit_BC12: Weibull fit
% ax1-3:    Handle to subplot axes 
function plot_fits_sep(xdata, Aip, C, fit_FN92, fit_BH05, fit_BC12, ax1, ax2, ax3, col)

x = 1:10:C;               % X for PL and WB


% Power law
BH = 10^fit_BH05(1).*x.^(fit_BH05(2));
semilogy(ax2, x, BH, 'Color', col);

% Weibull
BC = fit_BC12(1).*(x/fit_BC12(2)).^(fit_BC12(3)-2) .* exp(-((x/fit_BC12(2)).^fit_BC12(3)));
semilogy(ax3, x, BC, 'Color', col);

% Exponential

if length(Aip) == 1 && Aip == 0
    FN = exp(fit_FN92(1)).*exp(fit_FN92(2).*x);
    semilogy(ax1, x, FN, 'Color', col);
% 2 segments
elseif length(Aip) == 1 && Aip > 0
    [~, idx] = min(abs(x-sqrt(xdata(Aip))));
    x1 = x(1:idx);
    x2 = x(idx:end);
    
    FN1 = exp(fit_FN92(1,1)).*exp(fit_FN92(1,2).*x1);
    FN2 = exp(fit_FN92(2,1)).*exp(fit_FN92(2,2).*x2);
    
    semilogy(ax1, x1, FN1, 'Color', col);
    semilogy(ax1, x2, FN2, 'Color', col);
% 3 segments
else
    [~, idx1] = min(abs(x-sqrt(xdata(Aip(1)))));
    [~, idx2] = min(abs(x-sqrt(xdata(Aip(2)))));
    
    x1 = x(1:idx1);
    x2 = x(idx1:idx2);
    x3 = x(idx2:end);
    
    FN1 = exp(fit_FN92(1,1)).*exp(fit_FN92(1,2).*x1);
    FN2 = exp(fit_FN92(2,1)).*exp(fit_FN92(2,2).*x2);
    FN3 = exp(fit_FN92(3,1)).*exp(fit_FN92(3,2).*x3);
    
    semilogy(ax1, x1, FN1, 'Color', col);
    semilogy(ax1, x2, FN2, 'Color', col);
    semilogy(ax1, x3, FN3, 'Color', col);
end
