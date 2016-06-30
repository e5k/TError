% TError package
% Returns a matrix of relative and absolute errors
% Dist:     Type of error distribution
%   1: Uniform
%   2: Gaussian
% sims:     Number of simulations, i.e. size of the vector
% val:      Reference value
% err:      Relative error (%)
function mat = rand_err(dist, sims, val, err)

mat = zeros(sims, 2);   % Column 1: Errors in %
                        % Column 2: Values with errors
                        
if dist == 1            % Uniform error distribution
    mat(:,1) = 2 .* err .* rand(sims, 1) - err;    
    mat(:,2) = val + val .* mat(:,1)./100;
elseif dist == 2        % Gaussian error distribution 
    mat(:,1) = rand_G(0, err/3, sims);    
    mat(:,2) = val + val .* mat(:,1) ./100;
end
