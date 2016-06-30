% TERROR_sensitivity
% By: Sebastien Biass, Gholamhossein Bagheri, William Aeberhard and
% Costanza Bonadonna
% University of Geneva
% Copyright (C) 2014
%
%
% Email contact: costanza.bonadonna@unige.ch, sebastien.biasse@unige.ch
%
% This program is free software; 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation. 
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

function TError_sensitivity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_nm  = 'example';   % Run name

vent_ht = 6000;     % Vent elevation (m)
trop_ht = 17000;    % Tropopause height  (m)

rangeE  = -40:5:40; % Error vector applied to each input parameter sequentially

% Plume height - Carey and Sparks (1986)
dw_v    = 20.6;     % Downwind (km)
cw_v    = 11.3;     % Crosswind (km)
dm_v    = 1.6;      % Diameter (cm)
cl_d    = 2500;     % Clast density (kgm-3)

% MER - Wilson and Walker (1987)
cstWW_v = 0.295;    % Constant

% MER - Mastin et al. (2009)
cstMa_v = 2;        % Constant
%cstMa_v = 140.88;     % DRE density (kg/m3)
% MER - Degruyter and Bonadonna (2012)
wind_v  = -1;       % Average wind speed below tropopause (m/s)
                    % Set -1 to propagate the wind speed obtained from
                    % Carey and Sparks (1986)

% Volume
fl      = 'example.txt';
                    % Main file for volume calculation in a tab-delimited 
                    % text file where:
                    % Row 1: Location of the breaks-in-slope (BIS) for the
                    % exponential method, max of 3 segments (2 BIS). Enter
                    % 0 to use only 1 segment
                    % Row 2:n :
                    % Col 1: Thickness (cm)
                    % Col 2: Thickness error (%)
                    % Col 3: Area (km2)
                    % Col 4: Area error (%)
                    
% Volume - Bonadonna and Houghton (2005)
C_v     = 300;      % Distal integration limit (km)

% Volume - Bonadonna and Costa (2012)
% Ranges of lambda and n for the volume calculation using the Weibull fit.
% If ranges are entered as empty matrices (e.g. []), the code will
% automatically use the ranges suggested by Table 2 of Bonadonna and Costa
% (2013, Bulletin of Volcanology) using an average value of volume based on
% the exponential and power-law methods to estimate the VEI.
lam_r   = [5,100];   % Range of lambda for optimization algorithm
n_r     = [.5,50];   % Range of n for optimization algorithm
% NOTE: These ranges were defined for the example provided here, change for
% your own data

% Duration
dep_d_v = 1000;     % Bulk deposit density (kg/m3)

% Add path to functions
addpath(genpath('dep/'));

% Maximum error to plot
max_err = 1500;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beginning of calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstart = now;
clc;
display('_________________________________________________________________');
display(sprintf('TError_sensitivity run %s started: %s', run_nm, datestr(tstart)));


% Create the output folders
display('- Creating output folders...')
if exist(['Output/', run_nm, '/Sensitivity'], 'dir')
    choice = questdlg('The output folder already exists. Overwrite?', ...
        '', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            rmdir(['Output/', run_nm, '/Sensitivity'], 's');
            mkdir(['Output/', run_nm, '/Sensitivity/']);
        case 'No'
            display('Enter a different run name');
            return;
    end
else
    if exist(['Output/', run_nm, '/'], 'dir')
        mkdir(['Output/', run_nm, '/Sensitivity/']);
    else
        mkdir(['Output/', run_nm, '/']);                % Output folder
        mkdir(['Output/', run_nm, '/Sensitivity/']);    % Figure folder
    end
end


%% Volume parameters
file    = dlmread(fl);      

bis     = file(1,:); % Index value of the break in slope

file    = file(2:end, 1:4);
file    = flipud(sortrows(file, 1)); 
                     % Sort thickness in descending order

xdata   = file(:,3); % Area (km2)
ydata   = file(:,1); % Thickness (cm)

% Defines the breaks-in-slopes
if bis(1) == 0
    Aip = 0;
elseif bis(1) > 0 && bis(2) == 0
    Aip = bis(1);
elseif bis(2) > 0 && bis(3) == 0
    Aip = bis(1:2);
end

% Vector containing the different input parameters         
param = {'dw', 'cw', 'dm', 'const', 'dre', 'dep', 'wind', 'c', 'thick', 'area'};

% Storage matrix
stor = zeros(20, 10, length(rangeE),2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation using reference values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Running calculations...');

% Plume height and wind speed
[ht, wd]= get_height_CS86(dw_v, cw_v, dm_v, cl_d);      % Plume height, wind

% MER
mer(1)  = get_MER_WW87(ht - vent_ht/1000, cstWW_v);     % MER - Wilson & Walker
mer(2)  = get_MER_M09(ht - vent_ht/1000, cstMa_v);      % MER - Mastin et al.
if wind_v == -1
    mer(3) = get_MER_DB12(ht - vent_ht/1000, wd, (trop_ht-vent_ht)/1000);       % MER - Degruyter and Bonadonna with wind from Carey and Sparks
else
    mer(3) = get_MER_DB12(ht - vent_ht/1000, wind_v, (trop_ht-vent_ht)/1000);   % MER - Degruyter and Bonadonna with user-defined wind
end

% Volume
[vol(1), fit_FN92] = fn1992(xdata.^0.5, log(ydata), Aip);
[vol(2), ~] = bh2005(xdata.^0.5, ydata, fit_FN92(1), C_v);
[vol(3), ~] = bc2012(xdata.^0.5, ydata, lam_r, n_r);

% Mass
mass(1) = vol(1) * 10^9 * dep_d_v;
mass(2) = vol(2) * 10^9 * dep_d_v;
mass(3) = vol(3) * 10^9 * dep_d_v;

% Duration
dur(1)  = mass(1) / mer(1) / 60;
dur(2)  = mass(2) / mer(1) / 60;
dur(3)  = mass(3) / mer(1) / 60;
dur(4)  = mass(1) / mer(2) / 60;
dur(5)  = mass(2) / mer(2) / 60;
dur(6)  = mass(3) / mer(2) / 60;
dur(7)  = mass(1) / mer(3) / 60;
dur(8)  = mass(2) / mer(3) / 60;
dur(9)  = mass(3) / mer(3) / 60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loops applying the rangeE variable to input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(rangeE)
    display(sprintf('%4.0f %%', rangeE(i)));
    for j = 1:length(param)
        if j == 1
            dw = dw_v + dw_v * rangeE(i) / 100;
        else
            dw = dw_v;        
        end
        if j == 2
            cw = cw_v + cw_v * rangeE(i) / 100;
        else
            cw = cw_v;        
        end
        if j == 3
            dm = dm_v + dm_v * rangeE(i) / 100;
        else
            dm = dm_v;      
        end
        if j == 4
            cstWW = cstWW_v + cstWW_v * rangeE(i) / 100;
        else
            cstWW = cstWW_v;       
        end
        if j == 5
            cstMa = cstMa_v + cstMa_v * rangeE(i) / 100;
        else
            cstMa = cstMa_v;     
        end
        if j == 6
            dep_d = dep_d_v + dep_d_v * rangeE(i) / 100;
        else
            dep_d = dep_d_v;        
        end
        if j == 7
            wind = wind_v + wind_v * rangeE(i) / 100;
        else
            wind = wind_v;        
        end
        if j == 8
            C = C_v + C_v * rangeE(i) / 100;
        else
            C = C_v;    
        end
        if j == 9
            yerr = ydata + ydata .* rangeE(i) ./ 100;
        else
            yerr = ydata;      
        end
        if j == 10
            xerr = xdata + xdata .* rangeE(i) ./ 100;
        else
            xerr = xdata;   
        end
        
        % Plume height and wind speed
        [stor(1,j,i,1), stor(2,j,i,1)]= get_height_CS86(dw, cw, dm, cl_d);            % Plume height, wind
                
        % MER
        stor(3,j,i,1)  = get_MER_WW87(stor(1,j,i,1), cstWW);                    % MER - Wilson & Walker
        stor(4,j,i,1)  = get_MER_M09(stor(1,j,i,1), cstMa);                     % MER - Mastin et al.
        if wind_v == -1
            stor(5,j,i,1) = get_MER_DB12(stor(1,j,i,1) - vent_ht/1000, stor(2,j,i,1), (trop_ht-vent_ht)/1000);       % MER - Degruyter and Bonadonna with wind from Carey and Sparks
        else
            stor(5,j,i,1) = get_MER_DB12(stor(1,j,i,1) - vent_ht/1000, wind, (trop_ht-vent_ht)/1000);   % MER - Degruyter and Bonadonna with user-defined wind
        end

        % Volume
        [stor(6,j,i,1), fit_FN92] = fn1992(xerr.^0.5, log(yerr), Aip);
        [stor(7,j,i,1), ~] = bh2005(xerr.^0.5, yerr, fit_FN92(1), C);
        if isempty(lam_r)
            [lam_r, n_r]   = get_WBL(mean([vol(1), vol(2)]));
        end
        [stor(8,j,i,1), ~] = bc2012(xerr.^0.5, yerr, lam_r, n_r);

        % Mass
        stor(9,j,i,1) = stor(6,j,i,1) * 10^9 * dep_d;
        stor(10,j,i,1) = stor(7,j,i,1) * 10^9 * dep_d;
        stor(11,j,i,1) = stor(8,j,i,1) * 10^9 * dep_d;

        % Duration
        stor(12,j,i,1)  = stor(9,j,i,1) / stor(3,j,i,1) / 60;
        stor(13,j,i,1)  = stor(10,j,i,1) / stor(3,j,i,1) / 60;
        stor(14,j,i,1)  = stor(11,j,i,1) / stor(3,j,i,1) / 60;
        stor(15,j,i,1)  = stor(9,j,i,1) / stor(4,j,i,1) / 60;
        stor(16,j,i,1)  = stor(10,j,i,1) / stor(4,j,i,1) / 60;
        stor(17,j,i,1)  = stor(11,j,i,1) / stor(4,j,i,1) / 60;
        stor(18,j,i,1)  = stor(9,j,i,1) / stor(5,j,i,1) / 60;
        stor(19,j,i,1)  = stor(10,j,i,1) / stor(5,j,i,1) / 60;
        stor(20,j,i,1)  = stor(11,j,i,1) / stor(5,j,i,1) / 60;

        % Normalization
        stor(1,j,i,2) = (stor(1,j,i,1)-ht) / ht*100;
        stor(2,j,i,2) = (stor(2,j,i,1)-wd) / wd*100;
        stor(3,j,i,2) = (stor(3,j,i,1)-mer(1)) / mer(1)*100;
        stor(4,j,i,2) = (stor(4,j,i,1)-mer(2)) / mer(2)*100;
        stor(5,j,i,2) = (stor(5,j,i,1)-mer(3)) / mer(3)*100;
        stor(6,j,i,2) = (stor(6,j,i,1)-vol(1)) / vol(1)*100;
        stor(7,j,i,2) = (stor(7,j,i,1)-vol(2)) / vol(2)*100;
        stor(8,j,i,2) = (stor(8,j,i,1)-vol(3)) / vol(3)*100;
        stor(9,j,i,2) = (stor(9,j,i,1)-mass(1)) / mass(1)*100;
        stor(10,j,i,2) = (stor(10,j,i,1)-mass(2)) / mass(2)*100;
        stor(11,j,i,2) = (stor(11,j,i,1)-mass(3)) / mass(3)*100;
        stor(12,j,i,2) = (stor(12,j,i,1)-dur(1)) / dur(1)*100;
        stor(13,j,i,2) = (stor(13,j,i,1)-dur(2)) / dur(2)*100;
        stor(14,j,i,2) = (stor(14,j,i,1)-dur(3)) / dur(3)*100;
        stor(15,j,i,2) = (stor(15,j,i,1)-dur(4)) / dur(4)*100;
        stor(16,j,i,2) = (stor(16,j,i,1)-dur(5)) / dur(5)*100;
        stor(17,j,i,2) = (stor(17,j,i,1)-dur(6)) / dur(6)*100;
        stor(18,j,i,2) = (stor(18,j,i,1)-dur(7)) / dur(7)*100;
        stor(19,j,i,2) = (stor(19,j,i,1)-dur(8)) / dur(8)*100;
        stor(20,j,i,2) = (stor(20,j,i,1)-dur(9)) / dur(9)*100;
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preparing the outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

leg     = {'Downwind range', 'Crosswind range', 'Clast diameter', 'WW87 Constant', 'Ma09 Constant', 'Deposit density', 'Wind speed', 'Distal integration limit', 'Thickness measurement', 'Area of isopach contour'};

mod     = {'Height (CS86)', 'Wind (CS86)',...
            'MER (WW87)', 'MER (Ma09)', 'MER (DB12)',...
            'Volume (FN92)', 'Volume (BH05)', 'Volume (BC12)',...
            'Mass (FN92)', 'Mass (BH05)', 'Mass (BC12)',...
            'Duration (WW87-FN92)', 'Duration (WW87-BH05)','Duration (WW87-BC12)',...
            'Duration (Ma09-FN92)', 'Duration (Ma09-BH05)','Duration (Ma09-BC12)',...
            'Duration (DB12-FN92)', 'Duration (DB12-BH05)','Duration (DB12-BC12)'};

% Plot figure
display('Preparing figures...');
cmap = linspecer(10);
count = 1;
for i = 1:20
    tmp = shiftdim(squeeze(stor(i,:,:,2)),1);
    figure('Position', [465,431,655,252], 'Color', 'w'); hold on
    title(mod{i});
    xlabel('RIU (%)');
    ylabel('ROD (%)');
    for j = 1:10
        plot(rangeE, tmp(:,j), 'Color', cmap(j,:));
    end
    set(gca, 'Box', 'on', 'XTick', min(rangeE):10:max(rangeE));
    lim = get(gca, 'YLim');
    if lim(2) > max_err
        lim = [min(min(tmp)) max_err];
        set(gca, 'YLim', lim);
    end
    legend(leg, 'Box', 'on', 'Location', 'NorthEastOutside');
    saveas(gca, ['Output/', run_nm, '/Sensitivity/', num2str(count, '%02.0f'), '.eps']);
    %export_fig(gca, ['Output/', run_nm, '/Sensitivity/', num2str(count, '%02.0f'), '.eps']);
    close(gcf);
    count = count+1;
end

% Writing the report
display('Preparing the report...');
out_file = ['Output/', run_nm, '/sensitivity.xls'];
for i = 1:length(rangeE)
    sheet_name  = num2str(rangeE(i), '%2.0f');
    data_cell   = num2cell(stor(:,:,i,2));
    data_cell   = cellfun(@(x) num2str(x, '%4.0f'), data_cell, 'UniformOutput',0);
    out_mat     = [' ', leg; mod', data_cell];
    xlswrite(out_file, out_mat, sheet_name);
end

display(sprintf('TError_sensitivity run %s finished: %s (time elapsed: %3.0f min)', run_nm, datestr(now), etime(datevec(now),datevec(tstart))/60));
display('_________________________________________________________________');



