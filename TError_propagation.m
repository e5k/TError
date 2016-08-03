% TERROR_propagation
% By: Sebastien Biass, Gholamhossein Bagheri, William Aeberhard and
% Costanza Bonadonna
% University of Geneva
% Copyright (C) 2014
%
% Updates:  April 2016: Minor bug on plottig function
%
% Email contact: costanza.bonadonna@unige.ch, sebastien.biasse@unige.ch
%
% This program is free software; 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation. 
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

function TError_propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variable definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_nm  = 'example2';% Run name

vent_ht = 6000;     % Vent elevation (m)
trop_ht = 17000;    % Tropopause height  (m)

nb_sims = 10;     % Number of iteration of Monte Carlo simulation

error_d = 2;        % Error distribution
                    % 1: Uniform  -> error for each parameter is the maximum error
                    % 2: Gaussian -> error for each parameter is the 3*sigma (i.e. ~99% confidence interval)

% Plume height - Carey and Sparks (1986)
dw_v    = 20.6;     % Downwind (km)
dw_e    = 20;       % Downwind error (%)
cw_v    = 11.3;     % Crosswind (km)
cw_e    = 20;       % Crosswind error (%)
dm_v    = 1.6;      % Diameter (cm)
dm_e    = 20;       % Diameter error (%)

cl_d    = 2500;     % Clast density (kgm-3)

% MER - Wilson and Walker (1987)
cstWW_v = 0.295;    % Constant
cstWW_e = 20;       % Constant error (%)

% MER - Mastin et al. (2009)
cstMa_v = 2;        % Constant
cstMa_e = 20;       % Constant error (%)

% MER - Degruyter and Bonadonna (2012)
wind_v  = -1;     % Maximum wind speed below tropopause (m/s)
                    % Set -1 to propagate the wind speed obtained from
                    % Carey and Sparks (1986)
wind_e  = 20;       % Wind speed error (%)

% Volume
fl      = 'isopach_example.txt';
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
C_e     = 20;       % Distal integration limit  error (%)

% Volume - Bonadonna and Costa (2012)
% Ranges of lambda and n for the volume calculation using the Weibull fit.
% If ranges are entered as empty matrices (e.g. []), the code will
% automatically use the ranges suggested by Table 2 of Bonadonna and Costa
% (2013, Bull Volc) using an average value of volume based on
% the exponential and power-law methods to estimate the VEI.
lam_r   = [5,100];   % Range of lambda for optimization algorithm
n_r     = [.5,50];   % Range of n for optimization algorithm
% NOTE: These ranges were defined for the example provided here, change for
% your own data

% Duration
dep_d_v = 1000;     % Bulk deposit density (kg/m3)
dep_d_e = 20;       % Bulk deposit density error (%)

% Plotting options
plt     = 1;        % Plot figure? 0/1
frmt    = '.eps';   % Output image format
max_err = 500;      % Error limit (%) to plot above which all values will be comprised in the most extreme bin.

% Report options
pcile   = [2, 5, 95, 98];
                    % Percentile used in the report. Note that pairs 1-4
                    % and 2-3 should be symmetrical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beginning of calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstart = now;
clc;
display('_________________________________________________________________');
display(sprintf('TError run %s started: %s', run_nm, datestr(tstart)));


% Create the output folders
display('- Creating output folders...')
if exist(['Output/', run_nm, '/Propagation'], 'dir')
    choice = questdlg('The output folder already exists. Overwrite?', ...
        '', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            rmdir(['Output/', run_nm, '/Propagation'], 's');
            mkdir(['Output/', run_nm, '/Propagation/']);
        case 'No'
            display('Enter a different run name');
            return;
    end
else
    if exist(['Output/', run_nm, '/'], 'dir')
        mkdir(['Output/', run_nm, '/Propagation/']);
    else
        mkdir(['Output/', run_nm, '/']);                % Output folder
        mkdir(['Output/', run_nm, '/Propagation/']);    % Figure folder
    end
end

% Add path to functions
addpath(genpath('dep/'));   % Add path to required functions

%% Volume parameters
display('- Loading volume file...')
file    = dlmread(fl);      % Read volume file      
                            
bis     = file(1,:); % Index value of the break in slope

file    = file(2:end, 1:4); % Get data out of fle
file    = flipud(sortrows(file, 1)); 
                     % Sort thickness in descending order
                     
xdata   = file(:,3); % Area (km2)
ydata   = file(:,1); % Thickness (cm)

xerr    = file(:,4); % Area error (%)
yerr    = file(:,2); % Thickness error (%)

% Defines the breaks-in-slopes
if bis(1) == 0                      % Case 1 segment
    Aip = 0;
elseif bis(1) > 0 && bis(2) == 0    % Case 2 segments
    Aip = bis(1);
elseif bis(2) > 0 && bis(3) == 0    % Case 3 segments
    Aip = bis(1:2);
end

%% Create random values for each input parameter
display('- Creating randomness...')
dw_d    = rand_err(error_d, nb_sims, dw_v, dw_e);       % Downind range
cw_d    = rand_err(error_d, nb_sims, cw_v, cw_e);       % Crosswind range
dm_d    = rand_err(error_d, nb_sims, dm_v, dm_e);       % Clast diameter
cstWW_d = rand_err(error_d, nb_sims, cstWW_v, cstWW_e); % WW87 Constant
cstMa_d = rand_err(error_d, nb_sims, cstMa_v, cstMa_e); % Ma09 Constant
wind_d  = rand_err(error_d, nb_sims, wind_v, wind_e);   % Wind speed
C_d     = rand_err(error_d, nb_sims, C_v, C_e);         % Distal integration limit
dep_d_d = rand_err(error_d, nb_sims, dep_d_v, dep_d_e); % Bulk deposit density

% Thickness and area
xdata_d = zeros(size(xdata,1), nb_sims, 2);     
ydata_d = zeros(size(ydata,1), nb_sims, 2);

for i = 1:nb_sims
    if error_d == 1
        xdata_d(:, i, 1) = -xerr + (xerr - (-xerr)) .* rand(size(xerr,1), 1);
        ydata_d(:, i, 1) = -yerr + (yerr - (-yerr)) .* rand(size(yerr,1), 1);        
    else
        xdata_d(:, i, 1) = rand_G(0, xerr./2, size(xerr,1));   
        ydata_d(:, i, 1) = rand_G(0, yerr./2, size(yerr,1));  
    end
    xdata_d(:, i, 2) = xdata + xdata .* xdata_d(:, i, 1)./100;
    ydata_d(:, i, 2) = ydata + ydata .* ydata_d(:, i, 1)./100;
end


%% Storage matrices
% Col 1: Errors (%)
% Col 2: Values ± errors
display('- Creating storage matrices...')

% Plume height and wind speed
HT      = zeros(nb_sims, 2, 2);
                    % Dim 3:
                    % 1: Height
                    % 2: Wind speed
                    
% MER                    
MER     = zeros(nb_sims, 2, 3);
                    % Dim 3:
                    % 1: Wilson and Walker (1987)
                    % 2: Mastin et al. (2009)
                    % 3: Degtuyter and Bonadonna (2012)
                    
% Volume                    
VOL     = zeros(nb_sims, 2, 3);
                    % Dim 3:
                    % 1: Exponential
                    % 2: Power Law
                    % 3: Weibull
                    
% Fits                    
if bis(1) == 0
    FIT_EXP = zeros(nb_sims, 2, 1, 2);
elseif bis(1) > 0 && bis(2) == 0
    FIT_EXP = zeros(nb_sims, 2, 2, 2);
elseif bis(2) > 0 && bis(3) == 0
    FIT_EXP = zeros(nb_sims, 2, 3, 2);
end                 
                    % Dim 1: sims
                    % Dim 2:
                        % 1: T0
                        % 2: k
                    % Dim 3: segments
                    % Dim 4:
                        % 1: Values
                        % 2: Errors
FIT_PL  = zeros(nb_sims, 2, 2);
                    % Dim 3:
                    % 1: Values
                    % 2: Errors
FIT_WBL = zeros(nb_sims, 2, 3);                    
                    % Dim 3:
                    % 1: Values
                    % 2: Errors

% Mass
MASS    = zeros(nb_sims, 2, 3);
                    % Dim 3:
                    % 1: Exponential
                    % 2: Power Law
                    % 3: Weibull

% Duration
DUR     = zeros(nb_sims, 2, 9);
                    % Dim 3:
                    % 1: Wilson & Walker - Exponential
                    % 2: Wilson & Walker - Power Law
                    % 3: Wilson & Walker - Weibull
                    % 4: Mastin et al. - Exponential
                    % 5: Mastin et al. - Power Law
                    % 6: Mastin et al. - Weibull
                    % 7: Degruyter & Bonadonna - Exponential
                    % 8: Degruyter & Bonadonna - Power Law
                    % 9: Degruyter & Bonadonna - Weibull


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation using reference values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('- Calculations using reference values...')
% Plume height and wind speed
[ht, wd]= get_height_CS86(dw_v, cw_v, dm_v, cl_d);                          % Plume height, wind

% MER
mer(1)  = get_MER_WW87(ht - vent_ht/1000, cstWW_v);                         % MER - Wilson & Walker
mer(2)  = get_MER_M09(ht - vent_ht/1000, cstMa_v);                          % MER - Mastin et al.
if wind_v == -1
    mer(3) = get_MER_DB12(ht - vent_ht/1000, wd, (trop_ht-vent_ht)/1000);   % MER - Degruyter and Bonadonna with wind from Carey and Sparks
else
    mer(3) = get_MER_DB12(ht - vent_ht/1000, wind_v, (trop_ht-vent_ht)/1000);   % MER - Degruyter and Bonadonna with user-defined wind
end

% Volume
[vol(1), fit_FN92] = fn1992(xdata.^0.5, log(ydata), Aip);
[vol(2), fit_BH05] = bh2005(xdata.^0.5, ydata, fit_FN92(1), C_v);
if isempty(lam_r)
    [lam_r, n_r]   = get_WBL_ranges(mean([vol(1), vol(2)]));
end
[vol(3), fit_BC12] = bc2012(xdata.^0.5, ydata, lam_r, n_r);

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
%% Height, wind and MER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plume height and wind speed
display('- Propagating error in plume height...')
for i = 1:nb_sims
    [HT(i,2,1), HT(i,2,2)] = get_height_CS86(dw_d(i,2), cw_d(i,2), dm_d(i,2), cl_d);
end
HT(:,1,1)   = (HT(:,2,1) - ht) ./ ht*100;           % Normalize over median values
HT(:,1,2)   = (HT(:,2,2) - wd) ./ wd*100;           % Normalize over median values

% MER
display('- Propagating error in MER...')
MER(:,2,1)  = get_MER_WW87(HT(:,2,1) - vent_ht/1000, cstWW_d(:,2)); 
MER(:,2,2)  = get_MER_M09(HT(:,2,1) - vent_ht/1000, cstMa_d(:,2));
for i = 1:nb_sims
    if wind_v == -1
        MER(i,2,3) = get_MER_DB12(HT(i,2,1) - vent_ht/1000, HT(i,2,2), (trop_ht-vent_ht)/1000);
    else
        MER(i,2,3) = get_MER_DB12(HT(i,2,1) - vent_ht/1000, wind_d(i,2), (trop_ht-vent_ht)/1000); 
    end
end

MER(:,1,1) = (MER(:,2,1) - mer(1)) ./ mer(1)*100;   % Normalize over median values
MER(:,1,2) = (MER(:,2,2) - mer(2)) ./ mer(2)*100;   % Normalize over median values
MER(:,1,3) = (MER(:,2,3) - mer(3)) ./ mer(3)*100;   % Normalize over median values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume and Mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('- Propagating error in volume (this can take some time)...')
cmap = jet(nb_sims);
h2 = figure('Position', [560,340,562,608], 'Visible', 'off'); 
ax1 = subplot(3,1,1);
set(ax1, 'YScale', 'log');
title('Exponential');
xlabel('Area^0^.^5');
ylabel('Thickness (cm)');
xlim([0 C_v]);
ylim([10^-2 10^7]);
hold on

ax2 = subplot(3,1,2);
set(ax2, 'YScale', 'log');
title('Power Law');
xlabel('Area^0^.^5');
ylabel('Thickness (cm)');
xlim([0 C_v]);
ylim([10^-2 10^7]);
hold on 

ax3 = subplot(3,1,3);
set(ax3, 'YScale', 'log');
title('Weibull');
xlabel('Area^0^.^5');
ylabel('Thickness (cm)');
xlim([0 C_v]);
ylim([10^-2 10^7]);
hold on

% Fits volumes
wb = waitbar(0, 'Fitting isopach data...');
for i = 1:nb_sims
    [VOL(i,2,1), fit_tmp]    = fn1992(sqrt(xdata_d(:,i,2)), log(ydata_d(:,i,2)), Aip);    
    FIT_EXP(i,:,:,2) = reshape(fit_tmp', 1,2,size(fit_tmp,1)); 
    
    [VOL(i,2,2), fit_tmp_pl]     = bh2005(sqrt(xdata_d(:,i,2)), ydata_d(:,i,2), fit_tmp(1), C_d(i,2));
    FIT_PL(i,2,:) = reshape([10^fit_tmp_pl(1), fit_tmp_pl(2)],1,1,2);
    
    [VOL(i,2,3), fit_tmp_wbl]    = bc2012(sqrt(xdata_d(:,i,2)), ydata_d(:,i,2), lam_r, n_r);   
    FIT_WBL(i,2,:) = reshape(fit_tmp_wbl, 1,1,3);
    
    plot_fits_sep(xdata, Aip, C_d(i,2), fit_tmp, fit_tmp_pl, fit_tmp_wbl, ax1, ax2, ax3, cmap(i,:));
    
    waitbar(i/nb_sims);
end
close(wb);

hold off

VOL(:,1,1) = (VOL(:,2,1) - vol(1)) ./ vol(1)*100;
VOL(:,1,2) = (VOL(:,2,2) - vol(2)) ./ vol(2)*100;
VOL(:,1,3) = (VOL(:,2,3) - vol(3)) ./ vol(3)*100;

if bis(1) == 0
    FIT_EXP(:,1,1,1)  = (exp(FIT_EXP(:,1,1,2)) - exp(fit_FN92(1,1))) ./ exp(fit_FN92(1,1))*100;
    FIT_EXP(:,2,1,1)  = (FIT_EXP(:,2,1,2) - fit_FN92(1,2)) ./ fit_FN92(1,2)*100;
elseif bis(1) > 0 && bis(2) == 0
    FIT_EXP(:,1,1,1)  = (exp(FIT_EXP(:,1,1,2)) - exp(fit_FN92(1,1))) ./ exp(fit_FN92(1,1))*100;
    FIT_EXP(:,2,1,1)  = (FIT_EXP(:,2,1,2) - fit_FN92(1,2)) ./ fit_FN92(1,2)*100;
    FIT_EXP(:,1,2,1)  = (exp(FIT_EXP(:,1,2,2)) - exp(fit_FN92(2,1))) ./ exp(fit_FN92(2,1))*100;
    FIT_EXP(:,2,2,1)  = (FIT_EXP(:,2,2,2) - fit_FN92(2,2)) ./ fit_FN92(2,2)*100;
elseif bis(2) > 0 && bis(3) == 0
    FIT_EXP(:,1,1,1)  = (exp(FIT_EXP(:,1,1,2)) - exp(fit_FN92(1,1))) ./ exp(fit_FN92(1,1))*100;
    FIT_EXP(:,2,1,1)  = (fit_FN92(1,2) - FIT_EXP(:,2,1,2)) ./ fit_FN92(1,2)*100;
    FIT_EXP(:,1,2,1)  = (exp(FIT_EXP(:,1,2,2)) - exp(fit_FN92(2,1))) ./ exp(fit_FN92(2,1))*100;
    FIT_EXP(:,2,2,1)  = (fit_FN92(2,2) - FIT_EXP(:,2,2,2)) ./ fit_FN92(2,2)*100;
    FIT_EXP(:,1,3,1)  = (exp(FIT_EXP(:,1,3,2)) - exp(fit_FN92(3,1))) ./ exp(fit_FN92(3,1))*100;
    FIT_EXP(:,2,3,1)  = (fit_FN92(3,2) - FIT_EXP(:,2,3,2)) ./ fit_FN92(3,2)*100;
end

FIT_PL(:,1,1)   = (FIT_PL(:,2,1) - 10^fit_BH05(1)) ./ 10^fit_BH05(1)*100;
FIT_PL(:,1,2)   = (FIT_PL(:,2,2) - fit_BH05(2)) ./ fit_BH05(2)*100;

FIT_WBL(:,1,1)  = (FIT_WBL(:,2,1) - fit_BC12(1)) ./ fit_BC12(1)*100;
FIT_WBL(:,1,2)  = (FIT_WBL(:,2,2) - fit_BC12(2)) ./ fit_BC12(2)*100;
FIT_WBL(:,1,3)  = (FIT_WBL(:,2,3) - fit_BC12(3)) ./ fit_BC12(3)*100;

% Mass
display('- Propagating error in mass...')
MASS(:,2,1)     = VOL(:,2,1) .* 10^9 .* dep_d_d(:,2);
MASS(:,2,2)     = VOL(:,2,2) .* 10^9 .* dep_d_d(:,2);
MASS(:,2,3)     = VOL(:,2,3) .* 10^9 .* dep_d_d(:,2);

MASS(:,1,1)     = (MASS(:,2,1) - mass(1)) ./ mass(1)*100;
MASS(:,1,2)     = (MASS(:,2,2) - mass(2)) ./ mass(2)*100;
MASS(:,1,3)     = (MASS(:,2,3) - mass(3)) ./ mass(3)*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('- Propagating error in duration...')
DUR(:,2,1)      = MASS(:,2,1) ./ MER(:,2,1) ./ 60;
DUR(:,2,2)      = MASS(:,2,2) ./ MER(:,2,1) ./ 60;
DUR(:,2,3)      = MASS(:,2,3) ./ MER(:,2,1) ./ 60;
DUR(:,2,4)      = MASS(:,2,1) ./ MER(:,2,2) ./ 60;
DUR(:,2,5)      = MASS(:,2,2) ./ MER(:,2,2) ./ 60;
DUR(:,2,6)      = MASS(:,2,3) ./ MER(:,2,2) ./ 60;
DUR(:,2,7)      = MASS(:,2,1) ./ MER(:,2,3) ./ 60;
DUR(:,2,8)      = MASS(:,2,2) ./ MER(:,2,3) ./ 60;
DUR(:,2,9)      = MASS(:,2,3) ./ MER(:,2,3) ./ 60;

DUR(:,1,1)      = (DUR(:,2,1) - dur(1)) ./ dur(1)*100;
DUR(:,1,2)      = (DUR(:,2,2) - dur(2)) ./ dur(2)*100;
DUR(:,1,3)      = (DUR(:,2,3) - dur(3)) ./ dur(3)*100;
DUR(:,1,4)      = (DUR(:,2,4) - dur(4)) ./ dur(4)*100;
DUR(:,1,5)      = (DUR(:,2,5) - dur(5)) ./ dur(5)*100;
DUR(:,1,6)      = (DUR(:,2,6) - dur(6)) ./ dur(6)*100;
DUR(:,1,7)      = (DUR(:,2,7) - dur(7)) ./ dur(7)*100;
DUR(:,1,8)      = (DUR(:,2,8) - dur(8)) ./ dur(8)*100;
DUR(:,1,9)      = (DUR(:,2,9) - dur(9)) ./ dur(9)*100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Write the report
display('- Writing the output report...')
fid = fopen(['Output/', run_nm, '/', 'propagation.txt'], 'w');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
    'Parameter', 'Input value', 'Error', 'Mean', 'Median', 'Minimum', [num2str(pcile(1)),'th%'], [num2str(pcile(2)),'th%'], '25th%', '75th%', [num2str(pcile(3)),'th%'], [num2str(pcile(4)),'th%'], 'Maximum');


% Input parameters
writefile(fid, dw_d, dw_v, dw_e, 'Downwind range (km)', 0,0, pcile);
writefile(fid, cw_d, cw_v, cw_e, 'Crosswind range (km)', 0,0, pcile);
writefile(fid, dm_d, dm_v, dm_e, 'Clast diameter (cm)', 0,0, pcile);
writefile(fid, cstWW_d, cstWW_v, cstWW_e, 'WW87 constant', 0,0, pcile);
writefile(fid, cstMa_d, cstMa_v, cstMa_e, 'Ma09 constant', 0,0, pcile);
writefile(fid, dep_d_d, dep_d_v, dep_d_e, 'Deposit density (kgm-3)', 0,0, pcile);
if wind_v ~= -1
    writefile(fid, wind_d, wind_v, wind_e, 'Wind speed (ms-1)', 0,0, pcile);
end
writefile(fid, C_d, C_v, C_e, 'Distal integration limit (km)', 0,0, pcile);
fprintf(fid, '\n');

writefile(fid, HT(:,:,1), ht, [], 'Plume height (km asl)', 0,1, pcile);
writefile(fid, HT(:,:,2), wd, [], 'Wind speed (ms-1)', 0,1, pcile);
writefile(fid, MER(:,:,1), mer(1), [], 'MER - W&W87 (kgs-1)', 1,1, pcile);
writefile(fid, MER(:,:,2), mer(2), [], 'MER - Ma09 (kgs-1)', 1,1, pcile);
writefile(fid, MER(:,:,3), mer(3), [], 'MER - D&B12 (kgs-1)', 1,1, pcile);
fprintf(fid, '\n');

if bis(1) == 0
    writefile(fid, [FIT_EXP(:,1,1,1), exp(FIT_EXP(:,1,1,2))], exp(fit_FN92(1,1)), [], 'Exp-T01', 0,1, pcile);
    writefile(fid, [FIT_EXP(:,2,1,1), FIT_EXP(:,2,1,2)], fit_FN92(1,2), [], 'Exp-k1', 0,1, pcile);      
elseif bis(1) > 0 && bis(2) == 0
    writefile(fid, [FIT_EXP(:,1,1,1), exp(FIT_EXP(:,1,1,2))], exp(fit_FN92(1,1)), [], 'Exp-T01', 0,1, pcile);
    writefile(fid, [FIT_EXP(:,2,1,1), FIT_EXP(:,2,1,2)], fit_FN92(1,2), [], 'Exp-k1', 0,1, pcile);    
    writefile(fid, [FIT_EXP(:,1,2,1), exp(FIT_EXP(:,1,2,2))], exp(fit_FN92(2,1)), [], 'Exp-T02', 0,1, pcile);
    writefile(fid, [FIT_EXP(:,2,2,1), FIT_EXP(:,2,2,2)], fit_FN92(2,2), [], 'Exp-k2', 0,1, pcile);   
elseif bis(2) > 0 && bis(3) == 0
    writefile(fid, [FIT_EXP(:,1,1,1), exp(FIT_EXP(:,1,1,2))], exp(fit_FN92(1,1)), [], 'Exp-T01', 0,1, pcile);
    writefile(fid, [FIT_EXP(:,2,1,1), FIT_EXP(:,2,1,2)], fit_FN92(1,2), [], 'Exp-k1', 0,1, pcile);    
    writefile(fid, [FIT_EXP(:,1,2,1), exp(FIT_EXP(:,1,2,2))], exp(fit_FN92(2,1)), [], 'Exp-T02', 0,1, pcile);
    writefile(fid, [FIT_EXP(:,2,2,1), FIT_EXP(:,2,2,2)], fit_FN92(2,2), [], 'Exp-k2', 0,1, pcile);   
    writefile(fid, [FIT_EXP(:,1,3,1), exp(FIT_EXP(:,1,3,2))], exp(fit_FN92(3,1)), [], 'Exp-T03', 0,1, pcile);
    writefile(fid, [FIT_EXP(:,2,3,1), FIT_EXP(:,2,3,2)], fit_FN92(3,2), [], 'Exp-k3', 0,1, pcile);
end
writefile(fid, FIT_PL(:,:,1), 10^(fit_BH05(1)), [], 'PL-Tpl', 0,1, pcile);
writefile(fid, FIT_PL(:,:,2), fit_BH05(2), [], 'PL-k', 0,1, pcile);
writefile(fid, FIT_WBL(:,:,1), fit_BC12(1), [], 'WBL-Theta', 0,1, pcile);
writefile(fid, FIT_WBL(:,:,2), fit_BC12(2), [], 'WBL-Lambda', 0,1, pcile);
writefile(fid, FIT_WBL(:,:,3), fit_BC12(3), [], 'WBL-n', 0,1, pcile);
fprintf(fid, '\n');

writefile(fid, VOL(:,:,1), vol(1), [], 'Volume Exp (km3)', 1,1, pcile);
writefile(fid, VOL(:,:,2), vol(2), [], 'Volume PL (km3)', 1,1, pcile);
writefile(fid, VOL(:,:,3), vol(3), [], 'Volume WBL (km3)', 1,1, pcile);
fprintf(fid, '\n');

writefile(fid, MASS(:,:,1), mass(1), [], 'Mass Exp (km3)', 1,1, pcile);
writefile(fid, MASS(:,:,2), mass(2), [], 'Mass PL (km3)', 1,1, pcile);
writefile(fid, MASS(:,:,3), mass(3), [], 'Mass WBL (km3)', 1,1, pcile);
fprintf(fid, '\n');

writefile(fid, [DUR(:,1,1), DUR(:,2,1)], dur(1), [], 'Duration W&W87-Exp (h)', 0,1, pcile);
writefile(fid, [DUR(:,1,2), DUR(:,2,2)], dur(2), [], 'Duration W&W87-PL (h)', 0,1, pcile);
writefile(fid, [DUR(:,1,3), DUR(:,2,3)], dur(3), [], 'Duration W&W87-WBL (h)', 0,1, pcile);
writefile(fid, [DUR(:,1,4), DUR(:,2,4)], dur(4), [], 'Duration Ma09-Exp (h)', 0,1, pcile);
writefile(fid, [DUR(:,1,5), DUR(:,2,5)], dur(5), [], 'Duration Ma09-PL (h)', 0,1, pcile);
writefile(fid, [DUR(:,1,6), DUR(:,2,6)], dur(6), [], 'Duration Ma09-WBL (h)', 0,1, pcile);
writefile(fid, [DUR(:,1,7), DUR(:,2,7)], dur(7), [], 'Duration D&B12-Exp (h)', 0,1, pcile);
writefile(fid, [DUR(:,1,8), DUR(:,2,8)], dur(8), [], 'Duration D&B12-PL (h)', 0,1, pcile);
writefile(fid, [DUR(:,1,9), DUR(:,2,9)], dur(9), [], 'Duration D&B12-WBL (h)', 0,1, pcile);
fclose(fid);


%% Plot figures
if plt == 1
    display('- Plotting results and saving figures...')
    count = 1;
    % Input parameters
    plot_results(dw_d, dw_v, max_err, nb_sims, 'Downwind range (km)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_downwind', frmt]); count = count + 1;
    plot_results(cw_d, cw_v, max_err, nb_sims, 'Crosswind range (km)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_crosswind', frmt]); count = count + 1;
    plot_results(dm_d, dm_v, max_err, nb_sims, 'Clast diameter (cm)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_diameter', frmt]); count = count + 1;
    plot_results(cstWW_d, cstWW_v, max_err, nb_sims, 'WW87 Constant', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_WW87cst', frmt]); count = count + 1;
    plot_results(cstMa_d, cstMa_v, max_err, nb_sims, 'Ma09 Constant', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_Ma09cst', frmt]); count = count + 1;
    plot_results(dep_d_d, dep_d_v, max_err, nb_sims, 'Deposit density (kgm^-^3)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_deposit_density', frmt]); count = count + 1;
    plot_results(C_d, C_v, max_err, nb_sims, 'PL integration limit (km)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_C', frmt]); count = count + 1;
    if wind_v ~= -1
       plot_results(dw_d, dw_v, max_err, nb_sims, 'Wind speed (ms^-^1)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_wind_speed', frmt]); count = count + 1;
    end
    % Plume height
    plot_results(HT(:,:,1), ht, max_err, nb_sims, 'Plume height (km asl)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_plume_height_CS86', frmt]); count = count + 1;
    plot_results(HT(:,:,2), wd, max_err, nb_sims, 'Wind speed (ms^-^1)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_wind_speed_CS86', frmt]); count = count + 1;
    % MER
    plot_results(MER(:,:,1), mer(1), max_err, nb_sims, 'MER (kgs^-^1)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_MER_WW87', frmt]); count = count + 1;
    plot_results(MER(:,:,2), mer(2), max_err, nb_sims, 'MER (kgs^-^1)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_MER_Ma09', frmt]); count = count + 1;
    plot_results(MER(:,:,3), mer(3), max_err, nb_sims, 'MER (kgs^-^1)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_MER_DB12', frmt]); count = count + 1;
    % Volume fits
    plot_fits(xdata, ydata, xerr, yerr, Aip, C_v, fit_FN92, fit_BH05, fit_BC12); saveas(gcf, ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_fits', frmt]); close(gcf); count = count + 1;
    saveas(h2, ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_fits_detail', frmt]); close(h2); count = count + 1;
    if bis(1) == 0       
        plot_results([FIT_EXP(:,1,1,1), exp(FIT_EXP(:,1,1,2))], exp(fit_FN92(1,1)), max_err, nb_sims, 'T_01', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_T01_FN92', frmt]); count = count + 1;
        plot_results([FIT_EXP(:,2,1,1), FIT_EXP(:,2,1,2)], fit_FN92(1,2), max_err, nb_sims, 'k1', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_k1_FN92', frmt]); count = count + 1;      
    elseif bis(1) > 0 && bis(2) == 0
        plot_results([FIT_EXP(:,1,1,1), exp(FIT_EXP(:,1,1,2))], exp(fit_FN92(1,1)), max_err, nb_sims, 'T_01', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_T01_FN92', frmt]); count = count + 1;
        plot_results([FIT_EXP(:,2,1,1), FIT_EXP(:,2,1,2)], fit_FN92(1,2), max_err, nb_sims, 'k1', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_k1_FN92', frmt]); count = count + 1;    
        plot_results([FIT_EXP(:,1,2,1), exp(FIT_EXP(:,1,2,2))], exp(fit_FN92(2,1)), max_err, nb_sims, 'T_02', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_T02_FN92', frmt]); count = count + 1;
        plot_results([FIT_EXP(:,2,2,1), FIT_EXP(:,2,2,2)], fit_FN92(2,2), max_err, nb_sims, 'k2', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_k2_FN92', frmt]); count = count + 1;   
    elseif bis(2) > 0 && bis(3) == 0
        plot_results([FIT_EXP(:,1,1,1), exp(FIT_EXP(:,1,1,2))], exp(fit_FN92(1,1)), max_err, nb_sims, 'T_01', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_T01_FN92', frmt]); count = count + 1;
        plot_results([FIT_EXP(:,2,1,1), FIT_EXP(:,2,1,2)], fit_FN92(1,2), max_err, nb_sims, 'k1', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_k1_FN92', frmt]); count = count + 1;    
        plot_results([FIT_EXP(:,1,2,1), exp(FIT_EXP(:,1,2,2))], exp(fit_FN92(2,1)), max_err, nb_sims, 'T_02', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_T02_FN92', frmt]); count = count + 1;
        plot_results([FIT_EXP(:,2,2,1), FIT_EXP(:,2,2,2)], fit_FN92(2,2), max_err, nb_sims, 'k2', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_k2_FN92', frmt]); count = count + 1;   
        plot_results([FIT_EXP(:,1,3,1), exp(FIT_EXP(:,1,3,2))], exp(fit_FN92(3,1)), max_err, nb_sims, 'T_03', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_T03_FN92', frmt]); count = count + 1;
        plot_results([FIT_EXP(:,2,3,1), FIT_EXP(:,2,3,2)], fit_FN92(3,2), max_err, nb_sims, 'k3', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_k3_FN92', frmt]); count = count + 1;
    end
    plot_results(FIT_PL(:,:,1), 10^(fit_BH05(1)), max_err, nb_sims, 'T_p_l', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_Tpl_BH05', frmt]); count = count + 1;
    plot_results(FIT_PL(:,:,2), fit_BH05(2), max_err, nb_sims, 'k', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_k_BH05', frmt]); count = count + 1;
    plot_results(FIT_WBL(:,:,1), fit_BC12(1), max_err, nb_sims, '\theta', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_theta_BC12', frmt]); count = count + 1;
    plot_results(FIT_WBL(:,:,2), fit_BC12(2), max_err, nb_sims, '\lambda', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_lambda_BC12', frmt]); count = count + 1;
    plot_results(FIT_WBL(:,:,3), fit_BC12(3), max_err, nb_sims, 'n', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_n_BC12', frmt]); count = count + 1;
    % Volume
    plot_results(VOL(:,:,1), vol(1), max_err, nb_sims, 'Volume (km^3)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_volume_FN92', frmt]); count = count + 1;
    plot_results(VOL(:,:,2), vol(2), max_err, nb_sims, 'Volume (km^3)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_volume_BH05', frmt]); count = count + 1;
    plot_results(VOL(:,:,3), vol(3), max_err, nb_sims, 'Volume (km^3)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_volume_BC12', frmt]); count = count + 1;
    % Mass
    plot_results(MASS(:,:,1), mass(1), max_err, nb_sims, 'Mass (kg)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_mass_FN92', frmt]); count = count + 1;
    plot_results(MASS(:,:,2), mass(2), max_err, nb_sims, 'Mass (kg)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_mass_BH05', frmt]); count = count + 1;
    plot_results(MASS(:,:,3), mass(3), max_err, nb_sims, 'Mass (kg)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_mass_BC12', frmt]); count = count + 1;
    % Duration
    plot_results(DUR(:,:,1), dur(1), max_err, nb_sims, 'Duration (min)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_Duration_WW87-FN92', frmt]); count = count + 1;
    plot_results(DUR(:,:,2), dur(2), max_err, nb_sims, 'Duration (min)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_Duration_WW87-BH05', frmt]); count = count + 1;
    plot_results(DUR(:,:,3), dur(3), max_err, nb_sims, 'Duration (min)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_Duration_WW87-BC12', frmt]); count = count + 1;
    plot_results(DUR(:,:,4), dur(4), max_err, nb_sims, 'Duration (min)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_Duration_Ma09-FN92', frmt]); count = count + 1;
    plot_results(DUR(:,:,5), dur(5), max_err, nb_sims, 'Duration (min)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_Duration_Ma09-BH05', frmt]); count = count + 1;
    plot_results(DUR(:,:,6), dur(6), max_err, nb_sims, 'Duration (min)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_Duration_Ma09-BC12', frmt]); count = count + 1;
    plot_results(DUR(:,:,7), dur(7), max_err, nb_sims, 'Duration (min)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_Duration_DB12-FN92', frmt]); count = count + 1;
    plot_results(DUR(:,:,8), dur(8), max_err, nb_sims, 'Duration (min)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_Duration_DB12-BH05', frmt]); count = count + 1;
    plot_results(DUR(:,:,9), dur(9), max_err, nb_sims, 'Duration (min)', ['Output/', run_nm, '/Propagation/', num2str(count, '%02.0f'), '_Duration_DB12-BC12', frmt]);
end

display(sprintf('TError run %s finished: %s (time elapsed: %3.0f min)', run_nm, datestr(now), etime(datevec(now),datevec(tstart))/60));
display('_________________________________________________________________');


