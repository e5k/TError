% TError package
% This function writes the report of propagation runs
function writefile(fid, data, val_v, val_e, title, type1, type2, pcile)
% Type1 0: float
% Type1 1: power 10
data   = data(:,2);
tmp(1) = mean(data);
tmp(2) = prctile(data, 50); % Median
tmp(3) = min(data);         % Minimum
tmp(4) = prctile(data, pcile(1));  % 2th %
tmp(5) = prctile(data, pcile(2));  % 9th %
tmp(6) = prctile(data, 25); % 25th %

tmp(7) = prctile(data, 75); % 75th %
tmp(8) = prctile(data, pcile(3)); % 91th %
tmp(9) = prctile(data, pcile(4)); % 98th %
tmp(10) = max(data);         % Maximum



if type1 == 1
    fprintf(fid, '%s\t%6.1e\t%6.1e\t%6.1e\t%6.1e\t%6.1e\t%6.1e\t%6.1e\t%6.1e\t%6.1e\t%6.1e\t%6.1e\t%6.1e\n',...
        title, val_v, val_v*val_e/100, tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8), tmp(9), tmp(10));    
else
    fprintf(fid, '%s\t%4.1f\t%4.1f\t%4.1f\t%4.1f\t%4.1f\t%4.1f\t%4.1f\t%4.1f\t%4.1f\t%4.1f\t%4.1f\t%4.1f\n',...
        title, val_v, val_v*val_e/100, tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8), tmp(9), tmp(10));    
end

if type2 == 0
    fprintf(fid, '%s\t%s\t%4.1f%%\t%s\t%s\t%4.1f%%\t%4.1f%%\t%4.1f%%\t%4.1f%%\t%4.1f%%\t%4.1f%%\t%4.1f%%\t%4.1f%%\n',...
            '', '', val_e, '', '', (tmp(3)-tmp(2))/tmp(2)*100, (tmp(4)-tmp(2))/tmp(2)*100, (tmp(5)-tmp(2))/tmp(2)*100, (tmp(6)-tmp(2))/tmp(2)*100, (tmp(7)-tmp(2))/tmp(2)*100, (tmp(8)-tmp(2))/tmp(2)*100, (tmp(9)-tmp(2))/tmp(2)*100, (tmp(10)-tmp(2))/tmp(2)*100);
else
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%4.1f%%\t%4.1f%%\t%4.1f%%\t%4.1f%%\t%4.1f%%\t%4.1f%%\t%4.1f%%\t%4.1f%%\n',...
             '', '', '', '', '', (tmp(3)-tmp(2))/tmp(2)*100, (tmp(4)-tmp(2))/tmp(2)*100, (tmp(5)-tmp(2))/tmp(2)*100, (tmp(6)-tmp(2))/tmp(2)*100, (tmp(7)-tmp(2))/tmp(2)*100, (tmp(8)-tmp(2))/tmp(2)*100, (tmp(9)-tmp(2))/tmp(2)*100, (tmp(10)-tmp(2))/tmp(2)*100);
end