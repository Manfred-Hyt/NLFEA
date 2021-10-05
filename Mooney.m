%
% Calculate 2nd PK stress and stiffness for Mooney-Rivlin material
%
function [Stress D] = Mooney(F, A10, A01, K, ltan)
% Inputs:
%  F = Deformation gradient [3x3]
%  A10, A01, K = Material constants
%  ltan = 0 Calculate stress alone
%         1 Calculate stress and material stiffness
% Outputs:
%  Stress = 2nd PK stress [S11, S22, S33, S12, S23, S13];
%  D = Material stiffness [6x6]
%
X12 = 1/2; X13 = 1/3; X23 = 2/3; X43 = 4/3; X53 = 5/3; X89 = 8/9;
%
C = F'*F;
C1=C(1,1); C2=C(2,2); C3=C(3,3); C4=C(1,2); C5=C(2,3); C6=C(1,3);
I1 = C1+C2+C3;
I2 = C1*C2+C1*C3+C2*C3-C4^2-C5^2-C6^2;
I3 = det(C);
J3 = sqrt(I3);
J3M1 = J3 - 1.0D+00;
%
I1E = 2*[1 1 1 0 0 0]';
I2E = 2*[C2+C3, C3+C1, C1+C2, -C4, -C5, -C6]';
I3E = 2*[C2*C3-C5^2,  C3*C1-C6^2,  C1*C2-C4^2, ...
         C5*C6-C3*C4, C6*C4-C1*C5, C4*C5-C2*C6]';
%
W1 = I3^(-X13); W2 = X13*I1*I3^(-X43); W3 = I3^(-X23); 
W4 = X23*I2*I3^(-X53); W5 = X12*I3^(-X12);
%
J1E = W1*I1E - W2*I3E;
J2E = W3*I2E - W4*I3E;
J3E = W5*I3E;
%
Stress = A10*J1E + A01*J2E + K*J3M1*J3E;
%
% Material stiffness
%
D=zeros(6);
if ltan == 1
%
 I2EE = [0  4  4  0  0  0; 4  0  4  0  0  0; 4  4  0  0  0  0;
         0  0  0 -2  0  0; 0  0  0  0 -2  0; 0  0  0  0  0 -2];
 I3EE = [ 0     4*C3  4*C2  0    -4*C5  0;
          4*C3  0     4*C1  0     0    -4*C6;
          4*C2  4*C1  0    -4*C4  0     0;
          0     0    -4*C4 -2*C3  2*C6  2*C5;
         -4*C5  0     0     2*C6 -2*C1  2*C4;
          0    -4*C6  0     2*C5  2*C4 -2*C2];
%
 W1 = X23*I3^(-X12);    W2 = X89*I1*I3^(-X43); W3 = X13*I1*I3^(-X43);
 W4 = X43*I3^(-X12);    W5 = X89*I2*I3^(-X53); W6 = I3^(-X23);
 W7 = X23*I2*I3^(-X53); W8 = I3^(-X12);        W9 = X12*I3^(-X12);
%
 J1EE = -W1*(J1E*J3E' + J3E*J1E') + W2*(J3E*J3E') - W3*I3EE;
 J2EE = -W4*(J2E*J3E' + J3E*J2E') + W5*(J3E*J3E') + W6*I2EE - W7*I3EE;
 J3EE = -W8*(J3E*J3E') + W9*I3EE;
%
 D = A10*J1EE + A01*J2EE + K*(J3E*J3E') + K*J3M1*J3EE;
end
return;