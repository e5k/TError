% TError package
% Exponential fit for volume calulation with the method of Fierstein and Nathenson 1992
% xdata: Square root of the area (km)
% ydata: Log of thickness (cm)
% Aip:   Vector containing the break-in-slopes as indices
function [vol, fit_FN92_v]  = fn1992(xdata, ydata, Aip)

% 1 segment
if length(Aip) == 1 && Aip == 0
    N               = size(xdata,1);
    fit_FN92_v(2)   = ( sum(xdata.*ydata) - sum(xdata)*sum(ydata)/N ) / ( sum(xdata.^2) - sum(xdata)^2/N ); % Slope
    fit_FN92_v(1)   = mean(ydata)-fit_FN92_v(2)*mean(xdata);

% 2 segments
elseif length(Aip) == 1 && Aip > 0
    idx             = Aip;
    N               = size(xdata(1:idx),1);
    fit_FN92_v(1,2) = ( sum(xdata(1:idx).*ydata(1:idx)) - sum(xdata(1:idx))*sum(ydata(1:idx))/N ) / ( sum(xdata(1:idx).^2) - sum(xdata(1:idx))^2/N ); % Slope
    fit_FN92_v(1,1) = mean(ydata(1:idx))-fit_FN92_v(1,2)*mean(xdata(1:idx));
    
    N               = size(xdata(idx+1:end),1);
    fit_FN92_v(2,2) = ( sum(xdata(idx+1:end).*ydata(idx+1:end)) - sum(xdata(idx+1:end))*sum(ydata(idx+1:end))/N ) / ( sum(xdata(idx+1:end).^2) - sum(xdata(idx+1:end))^2/N ); % Slope
    fit_FN92_v(2,1) = mean(ydata(idx+1:end))-fit_FN92_v(2,2)*mean(xdata(idx+1:end));
 
% 3 segments
else
    idx1            = Aip(1);
    idx2            = Aip(2);   
    N               = size(xdata(1:idx1),1);
    fit_FN92_v(1,2) = ( sum(xdata(1:idx1).*ydata(1:idx1)) - sum(xdata(1:idx1))*sum(ydata(1:idx1))/N ) / ( sum(xdata(1:idx1).^2) - sum(xdata(1:idx1))^2/N ); % Slope
    fit_FN92_v(1,1) = mean(ydata(1:idx1))-fit_FN92_v(1,2)*mean(xdata(1:idx1));
      
    N               = size(xdata(idx1+1:idx2),1);
    fit_FN92_v(2,2) = ( sum(xdata(idx1+1:idx2).*ydata(idx1+1:idx2)) - sum(xdata(idx1+1:idx2))*sum(ydata(idx1+1:idx2))/N ) / ( sum(xdata(idx1+1:idx2).^2) - sum(xdata(idx1+1:idx2))^2/N ); % Slope
    fit_FN92_v(2,1) = mean(ydata(idx1+1:idx2))-fit_FN92_v(2,2)*mean(xdata(idx1+1:idx2)); 
    
    N               = size(xdata(idx2+1:end),1);
    fit_FN92_v(3,2) = ( sum(xdata(idx2+1:end).*ydata(idx2+1:end)) - sum(xdata(idx2+1:end))*sum(ydata(idx2+1:end))/N ) / ( sum(xdata(idx2+1:end).^2) - sum(xdata(idx2+1:end))^2/N ); % Slope
    fit_FN92_v(3,1) = mean(ydata(idx2+1:end))-fit_FN92_v(3,2)*mean(xdata(idx2+1:end));
    
end
 
vol                 = get_vol_FN92(exp(fit_FN92_v(:,1)),fit_FN92_v(:,2));