function [ X, T ] = Reaction_Sim( prop, x0, p, S, t_span )
%REACTION_SIM Stochastic simulation of chemical reactions
%   INPUT:
%   prop = function calculating propensities
%   x0 = row vector with starting conditions
%   p = row vector containing parameters
%   S = stoichiometry matrix
%   t = vector containing start and end time
%   OUTPUT:
%   T = row vector containing time steps
%   X = row vector containing amount of molecules for 
%   each species at time t
 
 t0 = t_span(1);
 t = t0;
 t_final = t_span(2);
 x = x0;
 
 counter = 1;
 n_var = length(x);
 X = zeros(10e6, n_var);
 X(counter, :) = x;
 T = zeros(1, 10e6);
 
while t < t_final
    u1 = rand();
    u2 = rand();
    propensities = prop(x, p);
    a0 = sum(propensities);
    tao = -log(u1)/a0;
    c = cumsum(propensities);
    r = find(c>=u2*a0, 1);
    nr = S(r,:);
    x = x + nr;
    t = t + tao;
    X(counter, :) = x;
    T(:, counter) = t;
    counter = counter + 1;
end
X = X(1:counter - 1, :);
T = T(:, 1:counter - 1);
end

