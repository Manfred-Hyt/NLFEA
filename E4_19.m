%
% Example 4.19 - Shear deformation of a square (finite rotation)
%
Young = 24000; nu=0.2; mu=Young/2/(1+nu); lambda=nu*Young/((1+nu)*(1-2*nu));
beta = 0; H = 1000; sY = 200*sqrt(3);
mp = [lambda mu beta H sY];
Iden=[1 1 1 0 0 0]';
D=2*mu*eye(6) + lambda*Iden*Iden';
D(4,4) = mu; D(5,5) = mu; D(6,6) = mu;
L = zeros(3,3);
stressN=[0 0 0 0 0 0]';
deps=[0 0 0 0 0 0]';
alphaN = [0 0 0 0 0 0]';
epN=0;
stressRN=stressN; alphaRN=alphaN;epRN=epN;
for i=1:15
    deps(4) = 0.004; L(1,2) = 0.024; L(2,1) = -0.02;
    [stressRN, alphaRN] = rotatedStress(L, stressRN, alphaRN);
    [stressR, alphaR, epR]=combHard(mp,D,deps,stressRN,alphaRN,epRN);
    [stress, alpha, ep]=combHard(mp,D,deps,stressN,alphaN,epN);
    X(i) = i*deps(4); Y1(i) = stress(4); Y2(i) = stressR(4);
    stressN = stress; alphaN = alpha; epN = ep;
    stressRN = stressR; alphaRN = alphaR; epRN = epR;
end
X = [0 X]; Y1=[0 Y1]; Y2=[0 Y2]; plot(X,Y1,X,Y2);