function w = prop_vilar(u, p)
%
% w = prop_vilar(u, p)
% Propensities, w, for the Vilar oscillator.
% 
% Input: u - the current state. 
%        p - list of parameters 
%
% The current state variables (u) are ordered as:
%    A  C  D_A  D_A'  D_R  D_R'  M_A  M_R  R
% The parameters (in p) are ordered as:
%   alfa_A  alfa'_A  alfa_R  alfa'_R  beta_A  beta_R  teta_A  teta_R ...
%   gamma_A  gamma_R  gamma_C  delta_M_R  delta_M_A  delta_A  delta_R
%

alfaA = p(1); alfapA = p(2);
alfaR = p(3);  alfapR = p(4);
betaA = p(5);  betaR = p(6);
tetaA = p(7);  tetaR = p(8);
gammaA = p(9); gammaR = p(10); gammaC = p(11);
deltaMR = p(12);  deltaMA = p(13); deltaA = p(14); 
deltaR = p(15);

w = zeros(18,1);
w(1 ) =   tetaA*u(4);
w(2 ) =  gammaA*u(1)*u(3);
w(3 ) =   tetaR*u(6);
w(4 ) =  gammaR*u(5)*u(1);
w(5 ) =  alfapR*u(6);
w(6 ) =   alfaR*u(5);
w(7 ) = deltaMR*u(8);
w(8 ) =  alfapA*u(4);
w(9 ) =   alfaA*u(3);
w(10) = deltaMA*u(7);
w(11) =   betaR*u(8);
w(12) =  deltaR*u(9);
w(13) =  deltaA*u(2);
w(14) =   betaA*u(7);
w(15) =   tetaA*u(4);
w(16) =   tetaR*u(6);
w(17) =  deltaA*u(1);
w(18) =  gammaC*u(1)*u(9);