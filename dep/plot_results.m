% TError package
% This function plots the results of the propagation run
function plot_results(data, in, maxERR, nb_sims, x_lab, out_name)

err     = data(:,1);
data    = data(:,2);
errmin  = min(err);             
errmax  = max(err);

% Sets plot boundaries with a maximum of maxERR
if errmax > maxERR
    maxX    = max(data(err<maxERR));
    errmax  = (maxX-in)/in*100;
else
    maxX    = max(data);
end

if errmin < -maxERR
    minX    = min(data(err>-maxERR));
    errmin  = (minX-in)/in*100;
else
    minX = min(data);
end

% Prepares the error vector to plot
if max([errmax, abs(errmin)])<10
    tmp = 2:2:10;    
elseif 10<max([errmax, abs(errmin)]) && max([errmax, abs(errmin)])<50
    tmp = 10:10:50;
elseif 50<max([errmax, abs(errmin)]) && max([errmax, abs(errmin)])<100
    tmp = 20:20:100;
elseif 100<max([errmax, abs(errmin)]) && max([errmax, abs(errmin)])<200
    tmp = 25:25:200;
elseif errmax>200
    tmp = 50:50:max([errmax, abs(errmin)])+50;
end
err_vec = [-fliplr(tmp), 0, tmp];

if errmax > abs(errmin)
   [~, idx] = min(abs(err_vec-errmin));
   if errmin < err_vec(idx)
       err_vec = err_vec(idx-1:end);
   else
       err_vec = err_vec(idx:end);
   end
else
   [~, idx] = min(abs(err_vec-errmax));
   if errmax > err_vec(idx)
       err_vec = err_vec(1:idx+1);
   else
       err_vec = err_vec(1:idx);
   end
    
end

if in < 0
    vec2plot = in + in.*-err_vec./100;
else
    vec2plot = in + in.*err_vec./100;
end

% Plot figure
figure('Position', [560,745,342,203], 'Color', 'w');
hold on;
for i = 1:length(err_vec)
    plot([vec2plot(i), vec2plot(i)], [0, 2*nb_sims], '-', 'Color', [.85 .85 .85], 'LineWidth', 1);
end
nhist(data, 'linewidth', .3, 'bplot', 'maxx', maxX, 'minx', minX, 'fsize', 10, 'ylabel', 'Frequency', 'xlabel', x_lab);
hax = gca;
axes('Position', get(hax, 'Position'), 'Color', 'none', 'XLim', get(hax, 'Xlim'), 'XAxisLocation','top', 'YAxisLocation','right', 'XTick', vec2plot, 'XTickLabel', strread(num2str(err_vec),'%s'), 'YTick', [], 'YTickLabel', []);
xlabel('ROD (%)')

saveas(gcf, out_name);
%export_fig(gcf, out_name);
close(gcf)
