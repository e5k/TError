% TError package
% In case the user did not specify it, this function returns typical ranges 
% of lambda and n values considering a VEI based on the mean of the volume 
% values obtained with the methods of Fierstein and Nathenson (1992) 
% and Bonadonna and Houghton (2005) to estimate the VEI.
% vol:  Volume (km3)
function [lam_r, n_r] = get_WBL_ranges(vol)
if vol >= 10
    lam_r   = [50, 1000];
    n_r     = [10, 100];
elseif vol < 10 && vol >= 1
    lam_r   = [20, 400];
    n_r     = [10, 100];
elseif vol < 1 && vol >= 0.1
    lam_r   = [5, 100];
    n_r     = [5, 50];
elseif vol < 0.1 && vol >= 0.01
    lam_r   = [1, 20];
    n_r     = [1, 10];
elseif vol < 0.01 && vol >= 0.001
    lam_r   = [0.5, 5];
    n_r     = [1, 10];
elseif vol < 0.001
    lam_r   = [0.2, 20];
    n_r     = [1, 10];
end
    
    