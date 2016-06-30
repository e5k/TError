% TError package
% Calculate the volume with the method of Bonadonna and Costa (2012)
% Theta, lambda, n: Fits obtained with the bc2012 function
function V = get_vol_BC12(theta, lambda, n)
% Theta:    Thickness scale (cm)
% Lambda:   Decay length scale of deposit thinning (km)
% n:        Shape parameter

% Equation 3 of Bonadonna and Costa (2012)
V = (2*theta/10^5*lambda^2)/n;       
