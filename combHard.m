%
% Linear combined isotropic/kinamtic hardening model
%
function [stress, alpha, ep]=combHard(mp,D,deps,stressN,alphaN,epN)
% Inputs:
% mp = [lambda, mu, beta, H, Y0];
% D = elastic stiffness matrix
% stressN = [s11, s22, s33, t12, t23, t13];
% alphaN = [a11, a22, a33, a12, a23, a13];
%
Iden = [1 1 1 0 0 0]'; 
two3 = 2/3; stwo3=sqrt(two3);                   %constants
mu=mp(2); beta=mp(3); H=mp(4); Y0=mp(5);        %material properties
ftol = Y0*1E-6;                                 %tolerance for yield
stresstr = stressN + D*deps;                    %trial stress
I1 = sum(stresstr(1:3));                        %trace(sigmatr)
str = stresstr - I1*Iden/3;                     %deviatoric stress
eta = str - alphaN;                             %shifted stress
etat = sqrt(eta(1)^2 + eta(2)^2 + eta(3)^2 ...
          + 2*(eta(4)^2 + eta(5)^2 + eta(6)^2));%norm of eta
fyld = etat - stwo3*(Y0+(1-beta)*H*epN);        %trial yield function
if fyld < ftol                                  %yield test
    stress = stresstr; alpha = alphaN; ep = epN;%trial states are final 
    return;
else
    gamma = fyld/(2*mu + two3*H);               %plastic consistency param
    ep = epN + gamma*stwo3;                     %updated eff. plastic strain
end
N = eta/etat;                                   %unit vector normal to f
stress = stresstr - 2*mu*gamma*N;               %updated stress
alpha = alphaN + two3*beta*H*gamma*N;           %updated back stress
