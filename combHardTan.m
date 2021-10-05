%
% Tangent stiffness for linear combined isotropic/kinamtic hardening model
%
function [Dtan]=combHardTan(mp,D,deps,stressN,alphaN,epN)
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
    Dtan = D; return;                           %elastic 
end
gamma = fyld/(2*mu + two3*H);                   %plastic consistency param
N = eta/etat;                                   %unit vector normal to f
var1 = 4*mu^2/(2*mu+two3*H);
var2 = 4*mu^2*gamma/etat;                       %coefficients
Dtan = D - (var1-var2)*N*N' + var2*Iden*Iden'/3;%tangent stiffness
Dtan(1,1) = Dtan(1,1) - var2;                   %contr. from 4th-order I
Dtan(2,2) = Dtan(2,2) - var2;
Dtan(3,3) = Dtan(3,3) - var2;
Dtan(4,4) = Dtan(4,4) - .5*var2;
Dtan(5,5) = Dtan(5,5) - .5*var2;
Dtan(6,6) = Dtan(6,6) - .5*var2;