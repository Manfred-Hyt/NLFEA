%
% 1D Linear combined isotropic/kinamtic hardening model
%
function [stress, alpha, ep]=combHard1D(mp,deps,stressN,alphaN,epN)
% Inputs:
% mp = [E, beta, H, Y0];
% deps = strain increment
% stressN = stress at load step N
% alphaN = back stress at load step N
% epN = plastic strain at load step N
%
E=mp(1); beta=mp(2); H=mp(3); Y0=mp(4);         %material properties
ftol = Y0*1E-6;                                 %tolerance for yield
stresstr = stressN + E*deps;                    %trial stress
etatr = stresstr - alphaN;                      %trial shifted stress
fyld = abs(etatr) - (Y0+(1-beta)*H*epN);        %trial yield function
if fyld < ftol                                  %yield test
    stress = stresstr; alpha = alphaN; ep = epN;%trial states are final 
    return;
else
    dep = fyld/(E+H);                           %plastic strain increment
end
stress = stresstr - sign(etatr)*E*dep;          %updated stress
alpha = alphaN + sign(etatr)*beta*H*dep;        %updated back stress
ep = epN + dep;                                 %updated plastic strain
return
