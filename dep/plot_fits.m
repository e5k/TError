% TError package
% This function plots the fits used for volume calculation obtained with
% reference values and adds a legend with the estimation of parameters
% xdata:    Square root of the area (km)
% ydata:    Thickness (cm)
% xerr:     Maximum error on square root of the area
% yerr:     Maximum error on thickness
% Aip:      Break-in-slopes as indices
% C:        Distal integration limit (km)
% fit_FN92: Exponential fit
% fit_BH05: Power-law fit
% fit_BC12: Weibull fit
function plot_fits(xdata, ydata, xerr, yerr, Aip, C, fit_FN92, fit_BH05, fit_BC12)

x = 1:10:C;               % X for PL and WB

figure;


% Power law
BH = 10^fit_BH05(1).*x.^(fit_BH05(2));
semilogy(x, BH, '-k', 'LineWidth', 2);
leg_BH = sprintf('Power Law:\t T_{PL} = %5.2f, k = %5.2f', 10^fit_BH05(1), fit_BH05(2));
hold on

% Weibull
BC = fit_BC12(1).*(x/fit_BC12(2)).^(fit_BC12(3)-2) .* exp(-((x/fit_BC12(2)).^fit_BC12(3)));
semilogy(x, BC, '-.k', 'LineWidth', 2);
leg_BC = sprintf('Weibull:\t \\theta = %5.2f, \\lambda = %5.2f, n = %5.2f', fit_BC12(1), fit_BC12(2), fit_BC12(3));

% Exponential

if length(Aip) == 1 && Aip == 0
    FN = exp(fit_FN92(1)).*exp(fit_FN92(2).*x);
    semilogy(x, FN, ':k', 'LineWidth', 2);
    leg_FN = sprintf('Exponential:\t T_{0} = %5.2f, k = %5.2f', exp(fit_FN92(1)), fit_FN92(2));
% 2 segments
elseif length(Aip) == 1 && Aip > 0
    [~, idx] = min(abs(x-sqrt(xdata(Aip))));
    x1 = x(1:idx);
    x2 = x(idx:end);
    
    FN1 = exp(fit_FN92(1,1)).*exp(fit_FN92(1,2).*x1);
    FN2 = exp(fit_FN92(2,1)).*exp(fit_FN92(2,2).*x2);
    
    leg_FN = sprintf('Exponential:\n\t 1: T_{01} = %5.2f, k_{1} = %5.2f\n\t 2: T_{02} = %5.2f, k_{2} = %5.2f', exp(fit_FN92(1,1)), fit_FN92(1,2), exp(fit_FN92(2,1)), fit_FN92(2,2));
    
    semilogy(x1, FN1, ':k', 'LineWidth', 2);
    semilogy(x2, FN2, ':k', 'LineWidth', 2);
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
    
    leg_FN = sprintf('Exponential:\n\t 1: T_{01} = %5.2f, k_{1} = %5.2f\n\t 2: T_{02} = %5.2f, k_{2} = %5.2f\n\t 3: T_{03} = %5.2f, k_{3} = %5.2f', exp(fit_FN92(1,1)), fit_FN92(1,2), exp(fit_FN92(2,1)), fit_FN92(2,2), exp(fit_FN92(3,1)), fit_FN92(3,2));
    
    semilogy(x1, FN1, ':k', 'LineWidth', 2);
    semilogy(x2, FN2, ':k', 'LineWidth', 2);
    semilogy(x3, FN3, ':k', 'LineWidth', 2);
end

xerr = sqrt(xdata + xdata .* xerr./100) - sqrt(xdata);
yerr = ydata .* yerr./100;

ploterr(sqrt(xdata), ydata, xerr, yerr, '.r', 'logy');

xlim([0 C]);
ylim([10^-2 10^ceil(max(log10(ydata)))]);
set(gca, 'Xlim', [0 C], 'YLim', [10^-2 10^(ceil(max(log10(ydata)))+1)]);
xlabel('Area^0^.^5');
ylabel('Thickness (cm)');
legend({leg_BH, leg_BC, leg_FN});

