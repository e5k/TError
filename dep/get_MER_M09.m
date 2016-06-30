% TError package
% Calculates the MER with the method of Mastin et al. (2009)
% ht:   Plume height (km above vent)
% const:Empirical constant
function MER = get_MER_M09(ht, const)

% Equation (1) of Mastin et al. (2009)
MER = ((ht./const).^(1/.241)).*2500;