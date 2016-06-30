% TError package
% Calculate plume height (km asl) and wind speed (ms) using an
% implementation of the model of Carey and Sparks (1986)
% dw:   Downwind range (km)
% cw:   Crosswind range (km)
% d:    Clast diameter (cm)
% den:  Clast density (kgm-3)
function [height, wind] = get_height_CS86(dw, cw, d, den)
% 1. Find values of height and wind for a given set of cw/dw on all Figures
%    16 of Carey and Sparks (1986)
% 2. Interpolates the values of height and wind for a given
%    diameter/density

d = d/100;  % Clast diameter in m
% 
% dw=20.;
% cw=10.;
% d=8e-2;
% den=2000.;

%% 1. Find values of height and wind for a given set of cw/dw on all Figures 16 of Carey and Sparks (1986)
% Figure 16a
d_a=0.8e-2;
den_a=2500;

a0=dw;
% Polynomial fit
a10=-9e-5*dw^3+0.0117*dw^2+.4887*dw;
a20=-8e-5*dw^3+0.0129*dw^2+.252*dw;
a30=-6e-5*dw^3+0.0114*dw^2+.1229*dw;

if (cw>=a0)
    wind_a=0;
elseif (cw<a0) && (cw>=a10)
    wind_a=0.+(10-0)*(cw-a0)/(a10-a0);
elseif (cw<a10) && (cw>=a20)
    wind_a=10.+(20-10)*(cw-a10)/(a20-a10);
elseif (cw<a20)
    wind_a=20.+(30-20)*(cw-a20)/(a30-a20); 
end
height_a=(-0.01977*cw^3 + 1.524*cw^2 + 16.26*cw + 0.06977) / (cw + 1.928);

% Figure 16b
d_b=1.6e-2;
den_b=2500;
b0=dw;
% Polynomial fit
b10=-2e-5*dw^3+0.0092*dw^2+.5168*dw;
b20=-10e-5*dw^3+0.0151*dw^2+.2577*dw;
b30=-10e-5*dw^3+0.0153*dw^2+.1191*dw;

if (cw>=b0)
    wind_b=0;
elseif (cw<b0) && (cw>=b10)
    wind_b=0.+(10-0)*(cw-b0)/(b10-b0);
elseif (cw<b10) && (cw>=b20)
    wind_b=10.+(20-10)*(cw-b10)/(b20-b10);
elseif (cw<b20)
    wind_b=20.+(30-20)*(cw-b20)/(b30-b20); 
end
height_b=(-0.02357*cw^3 + 1.672*cw^2 + 17.69*cw + 0.08577) / (cw + 1.557);   

% Figure 16c
d_c=3.2e-2;
den_c=2500;
c0=dw;
% Rational fit
c10=(0.01254 *dw^5 + 0.2165*dw^4 -9.467 *dw^3 +  85.12*dw^2 -199.4*dw -15.95) /(dw^3 -23.58*dw^2 +  183.3*dw -415);
c20= (0.01259*dw^5 +  0.0369 *dw^4 -5.405*dw^3 +55.24*dw^2 -104.7*dw + 22.13) /(dw^3 -23.97*dw^2 + 183.7*dw -294.2);
c30= (-0.0314*dw^5 + 13.41*dw^4  -85.39*dw^3  -4602*dw^2 + 5.544e+04*dw -3507 ) /(dw^3 + 1061*dw^2  -3.068e+04*dw + 2.413e+05);

if (cw>=c0)
    wind_c=0;
elseif (cw<c0) && (cw>=c10)
    wind_c=0.+(10-0)*(cw-c0)/(c10-c0);
elseif (cw<c10) && (cw>=c20)
    wind_c=10.+(20-10)*(cw-c10)/(c20-c10);
elseif (cw<c20)
    wind_c=20.+(30-20)*(cw-c20)/(c30-c20); 
end
height_c=(-0.03649*cw^3 +2.188*cw^2 + 15.27*cw -0.03349) / (cw +   0.8307);

% Figure 16d
d_d=6.4e-2;
den_d=2500;
d0=dw;
% Polynomial fit
d10=-0.0009*dw^3+0.0379*dw^2+.4446*dw;
d20=-0.0007*dw^3+0.04*dw^2+.1919*dw;
d30=-0.0006*dw^3+0.0354*dw^2+.0867*dw;

if (cw>=d0)
    wind_d=0;
elseif (cw<d0) && (cw>=d10)
    wind_d=0.+(10-0)*(cw-d0)/(d10-d0);
elseif (cw<d10) && (cw>=d20)
    wind_d=10.+(20-10)*(cw-d10)/(d20-d10);
elseif (cw<d20)
    wind_d=20.+(30-20)*(cw-d20)/(d30-d20); 
end
height_d=(-0.05547*cw^3 + 2.729*cw^2 + 18.71*cw + 0.001767) / (cw + 0.9346);


% 2. Interpolates the values of height and wind for a given diameter/density
if (d*den<=d_a*den_a)
    height=height_a+(height_b-height_a)*(d*den-d_a*den_a)/(d_b*den_b-d_a*den_a);
    wind=wind_a+(wind_b-wind_a)*(d*den-d_a*den_a)/(d_b*den_b-d_a*den_a);
elseif (d*den>=d_a*den_a)&& (d*den<d_b*den_b)
    height=height_a+(height_b-height_a)*(d*den-d_a*den_a)/(d_b*den_b-d_a*den_a);
    wind=wind_a+(wind_b-wind_a)*(d*den-d_a*den_a)/(d_b*den_b-d_a*den_a);
elseif (d*den>=d_b*den_b)&& (d*den<d_c*den_c)
    height=height_b+(height_c-height_b)*(d*den-d_b*den_b)/(d_c*den_c-d_b*den_b);
    wind=wind_b+(wind_c-wind_b)*(d*den-d_b*den_b)/(d_c*den_c-d_b*den_b);
elseif (d*den>=d_c*den_c)   
    height=height_c+(height_d-height_c)*(d*den-d_c*den_c)/(d_d*den_d-d_c*den_c);
    wind=wind_c+(wind_d-wind_c)*(d*den-d_c*den_c)/(d_d*den_d-d_c*den_c);
end
