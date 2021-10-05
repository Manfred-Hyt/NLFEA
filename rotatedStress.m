%
% Rotate stress and back stress to the mid-configuration
%
function [stress, alpha] = rotatedStress(L, S, A)
%L = [dui/dxj]
str=[S(1) S(4) S(6);S(4) S(2) S(5);S(6) S(5) S(3)];
alp=[A(1) A(4) A(6);A(4) A(2) A(5);A(6) A(5) A(3)];
R = L*inv(eye(3) - 0.5*L);
W = .5*(R-R');
R = eye(3) + inv(eye(3) - 0.5*W)*W;
str = R*str*R';
alp = R*alp*R';
stress=[str(1,1) str(2,2) str(3,3) str(1,2) str(2,3) str(1,3)]';
alpha =[alp(1,1) alp(2,2) alp(3,3) alp(1,2) alp(2,3) alp(1,3)]';
return;