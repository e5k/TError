% TError package
% Calculate the volume with the method of Fierstein and Nathenson (1992)
% T0    Extrapolated thickness at A = 0 (cm)
% k     Slope of the segment
function V = get_vol_FN92(T0, k)

T0 = T0/10^5;

%% 1 segment
if size(T0,1) == 1
    % Equation (12) of Fierstein and Nathenson (1992) 
    V = 2*(T0)/k^2;       
end    

%% 2 segments
if size(T0,1) == 2
    A = (log(T0(2))-log(T0(1)))/(k(1)-k(2));
    k = -k;
    % Equation (18) of Fierstein and Nathenson (1992) 
    V = 2*T0(1)/k(1)^2 + 2*T0(1)*((k(2)*A+1)/k(2)^2 - (k(1)*A+1)/k(1)^2)*exp(-k(1)*A);                % Equation 18 of Fierstein and Nathenson (1992)
end

%% 3 segments   
if size(T0,1) == 3
    A1 = (log(T0(2))-log(T0(1)))/(k(1)-k(2));
    A2 = (log(T0(3))-log(T0(2)))/(k(2)-k(3));
    k = -k;
    % Equation (3) of Bonadonna and Houghton (2005)
    V = 2*T0(1)/k(1)^2 + 2*T0(1)*((k(2)*A1+1)/k(2)^2 - (k(1)*A1+1)/k(1)^2)*exp(-k(1)*A1) + 2*T0(2)*((k(3)*A2+1)/k(3)^2 - (k(2)*A2+1)/k(2)^2)*exp(-k(2)*A2);
end
