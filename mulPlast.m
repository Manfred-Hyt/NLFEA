%
% Multiplicative plasticity with linear combined hardening
%
function [stress, b, alpha, ep]=mulPlast(mp,D,L,b,alpha,ep)
%mp = [lambda, mu, beta, H, Y0];
%D = elasticity matrix b/w prin stress & log prin stretch (3x3)
%L = [dui/dxj] velocity gradient
%b = elastic left C-G deformation vector (6x1)
%alpha = principal back stress (3x1)
%ep = effective plastic strain
%
EPS=1E-12;
Iden = [1 1 1]'; two3 = 2/3; stwo3=sqrt(two3); %constants
mu=mp(2); beta=mp(3); H=mp(4); Y0=mp(5);       %material properties
ftol = Y0*1E-6;                                %tolerance for yield
R = inv(eye(3)-L);                             %inc. deformation gradient
bm=[b(1) b(4) b(6);b(4) b(2) b(5);b(6) b(5) b(3)];
bm = R*bm*R';                                  %trial elastic left C-G
b=[bm(1,1) bm(2,2) bm(3,3) bm(1,2) bm(2,3) bm(1,3)]';
[~,P]=eig(bm);                                 %eigenvalues
eigen=sort(real([P(1,1) P(2,2) P(3,3)]))';           %principal stretch
%
% Duplicated eigenvalues
TMP=-1;
for I=1:2
  if abs(eigen(1)-eigen(3)) < EPS
    eigen(I)=eigen(I)+TMP*EPS;
    TMP=-TMP;
  end
end
if abs(eigen(1)-eigen(2)) < EPS; eigen(2) = eigen(2) + EPS; end;
if abs(eigen(2)-eigen(3)) < EPS; eigen(2) = eigen(2) + EPS; end;
%
% EIGENVECTOR MATRIX N*N' = M(6,*)
M=zeros(6,3);                                  %eigenvector matrices
for K=1:3
  KB=1+mod(K,3);
  KC=1+mod(KB,3);
  EA=eigen(K);
  EB=eigen(KB);
  EC=eigen(KC);
  D1=EB-EA;
  D2=EC-EA;
  DA = 1 / (D1 * D2);
  M(1,K)=((b(1)-EB)*(b(1)-EC)+b(4)*b(4)+b(6)*b(6))*DA;
  M(2,K)=((b(2)-EB)*(b(2)-EC)+b(4)*b(4)+b(5)*b(5))*DA;
  M(3,K)=((b(3)-EB)*(b(3)-EC)+b(5)*b(5)+b(6)*b(6))*DA;
  M(4,K)=(b(4)*(b(1)-EB+b(2)-EC)+b(5)*b(6))*DA;
  M(5,K)=(b(5)*(b(2)-EB+b(3)-EC)+b(4)*b(6))*DA;
  M(6,K)=(b(6)*(b(3)-EB+b(1)-EC)+b(4)*b(5))*DA;
end
%
eigen=sort(real([P(1,1) P(2,2) P(3,3)]))';           %principal stretch
deps = 0.5*log(eigen);                         %logarithmic 
sigtr = D*deps;                                %trial principal stress
eta = sigtr - alpha - sum(sigtr)*Iden/3;       %shifted stress
etat = norm(eta);                              %norm of eta
fyld = etat - stwo3*(Y0+(1-beta)*H*ep);        %trial yield function
if fyld < ftol                                 %yield test
   sig = sigtr;                                %trial states are final 
   stress = M*sig;                             %stress (6x1)
else
   gamma = fyld/(2*mu + two3*H);               %plastic consistency param
   ep = ep + gamma*stwo3;                      %updated eff. plastic strain
   N = eta/etat;                               %unit vector normal to f
   deps = deps - gamma*N;                      %updated elastic strain
   sig = sigtr - 2*mu*gamma*N;                 %updated stress
   alpha = alpha + two3*beta*H*gamma*N;        %updated back stress
   stress = M*sig;                             %stress (6x1)
   b = M*exp(2*deps);                          %updated elastic left C-G
end