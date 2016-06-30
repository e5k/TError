% TError package
% Calculates the MER with the method of Wilson and Walker (1987)
% ht:   Plume height (km above vent)
% const:Empirical constant
function MER = get_MER_WW87(Ht, cons)

% Equation (16) of Wilson and Walker (1987)
MER = (Ht./cons).^4;