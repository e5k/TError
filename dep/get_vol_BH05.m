% TError package
% Calculate the volume with the method of Bonadonna and Houghton (2005)
% Theta, lambda, n: Fits obtained with the bc2012 function
% T0:    Maximum thickness
% m:     Power law exponent
% TPl:   Coefficient (= Tpl)
% C:     Distal integration limit (km)
function V = get_vol_BH05(T0, m, Tpl, C)

T0    = T0/10^5;
Tpl   = Tpl/10^5;

% Equation 6 of Bonadonna and Houghton (2005)
V = (2*Tpl/(2-m)) * (C^(2-m) - ((T0/Tpl)^(-1/m))^(2-m));       
