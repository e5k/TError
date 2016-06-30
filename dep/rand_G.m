% TError package
% This function returns a Gausian noise
% mu:       Median
% sigma:    Standard deviation
% nb_runs:  Number of runs of Monte Carlo simulations
function r = rand_G(mu, sigma, nb_runs)
r = mu + sigma.*randn(nb_runs, 1);