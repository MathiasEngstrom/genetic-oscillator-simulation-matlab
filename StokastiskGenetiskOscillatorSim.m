% Script running a stochastic simulation of a circadian clock using the
% function Reaction_Sim

% Time span
t_span = [0 200];

% Setting parameters of the model and loading them into vector
alphaA = 50;
alphaAp = 500;
alphaR = 0.01;
alphaRp = 50;
betaA = 50;
betaR = 5;
deltaMA = 10;
deltaMR = 0.5;
deltaA = 1;
deltaR = 0.2; %0.2 resp 0.08
gammaA = 1;
gammaR = 1;
gammaC = 2;
thetaA = 50;
thetaR = 100;
p = [alphaA alphaAp alphaR alphaRp betaA betaR thetaA thetaR gammaA... 
    gammaR gammaC deltaMR deltaMA deltaA deltaR];

% Setting initial conditions and loading into vector
A = 0;
C = 0;
DA = 1;
DAp = 0;
DR = 1;
DRp = 0;
MA = 0;
MR = 0;
R = 0;
x0 = [A C DA DAp DR DRp MA MR R];

% Defining stoichiometry matrix
S = nr_vilar;

% Simulating the reaction, @prop_vilar is a function that calculates 
% the propensities
[ X, T ] = Reaction_Sim( @prop_vilar, x0, p, S, t_span );

% Plot amount of protein A and R against time
plot(T, [X(:, 1) X(:, 9)]);
  xlabel('Time [h]')
  ylabel('Number of protein molecules');
  title('Stochastic simulation of the circadian clock (deltaR = 0.2)');
  legend('A','R');