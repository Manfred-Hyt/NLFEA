%
% Example 4.22 - Shear deformation of a square
%
clear
Young = 24000; nu=0.2; mu=Young/2/(1+nu); lambda=nu*Young/((1+nu)*(1-2*nu));
beta = 0; H = 1000; sY = 200*sqrt(3);
mp = [lambda mu beta H sY];
Iden=[1 1 1 0 0 0]';
D=2*mu*eye(6) + lambda*Iden*Iden';
D(4,4) = mu; D(5,5) = mu; D(6,6) = mu;
Iden=[1 1 1]';
DM=2*mu*eye(3) + lambda*Iden*Iden';
L = zeros(3,3);
stressN=[0 0 0 0 0 0]';
deps=[0 0 0 0 0 0]';
alphaN = [0 0 0 0 0 0]';
epN=0;
stressRN=stressN; alphaRN=alphaN;epRN=epN;
bMN=[1 1 1 0 0 0]';
alphaMN = [0 0 0]';
epMN=0;
for i=1:15
    deps(4) = 0.004; L(1,2) = 0.024; L(2,1) = -0.02;
    [stressRN, alphaRN] = rotatedStress(L, stressRN, alphaRN);
    [stressR, alphaR, epR]=combHard(mp,D,deps,stressRN,alphaRN,epRN);
    [stress, alpha, ep]=combHard(mp,D,deps,stressN,alphaN,epN);
    [stressM, bM, alphaM, epM]=mulPlast(mp,DM,L,bMN,alphaMN,epMN);
    X(i)=i*deps(4);Y1(i)=stress(4);Y2(i)=stressR(4);Y3(i)=stressM(4);
    stressN = stress; alphaN = alpha; epN = ep;
    stressRN = stressR; alphaRN = alphaR; epRN = epR;
    bMN=bM; alphaMN = alphaM; epMN = epM;
end
X = [0 X]; Y1=[0 Y1]; Y2=[0 Y2]; Y3 = [0 Y3]; plot(X,Y1,X,Y2,X,Y3);