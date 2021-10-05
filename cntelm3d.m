function [FORCE, STIFF] = cntelm3d(OMEGAN, ELXY, LTAN)
%***********************************************************************
% SEARCH CONTACT POINT AND RETURN STIFFNESS AND RESIDUAL FORCE 
% IF CONTACTED FOR NORMAL CONTACT
%***********************************************************************
%
  EPS=1.E-6; TL1=2; TL2=0.1; TL3=1.01; FORCE=[]; STIFF=[];
%
% COMPUTE TWO TANGENT VECTORS & APPROXIMATED CONTACT POINT AT THE CENTER
  [T1, T2, XS] = CUTL(0,0,ELXY);
  XN = cross(T1, T2); XN=XN/norm(XN);
  XI  = (XS'*T1)/(2*norm(T1)^2);
  ETA = (XS'*T2)/(2*norm(T2)^2);
  GN = XN'*XS;
  XX=(ELXY(:,2)-ELXY(:,3)+ELXY(:,4)-ELXY(:,5))/4;
%
% IF NATURAL COORD. IS OUT OF BOUND, NO CONTACT IS ASSUMED
  if((XI<-TL1)||(XI>TL1)||(ETA<-TL1)||(ETA>TL1)||GN>TL2), return; end
%
% FIND EXACT CONTACT POINT THROUGH NEWTON-RAPHSON METHOD
  for ICOUNT=1:20
    [T1, T2, XS] = CUTL(ETA,XI,ELXY);
    A=[-T1'*T1, XS'*XX-T2'*T1; XS'*XX-T2'*T1, -T2'*T2];
    B=[-XS'*T1; -XS'*T2];
    DXI=A\B;
    XI=XI+DXI(1); ETA=ETA+DXI(2);
    if(norm(DXI)<EPS), break; end
  end
%
% CHECK FOR REFERENCE COORDINATE WITHIN RANGE
  if((XI<-TL3)||(XI>TL3)||(ETA<-TL3)||(ETA>TL3)), return; end
%
% NORMAL GAP FUNCTION
  XN = cross(T1, T2); XN=XN/norm(XN);
  GN = XN'*XS;
  if GN>0, return; end
%
% CONTACT FORCE
  FORCE = -OMEGAN*GN*XN;
%
% FORM STIFFNESS (NONFRICTION)
  if LTAN, STIFF = OMEGAN*(XN*XN'); end
end
function [T1, T2, XS] = CUTL(ETA,XI,ELXY)
%***********************************************************************
% COMPUTE COORD. OF CENTEROID AND TWO TANGENT VECTORS
%***********************************************************************
  XNODE=[0 -1 1 1 -1; 0 -1 -1 1 1];
  T1 = zeros(3,1); T2 = zeros(3,1); XC = zeros(3,1); XS = zeros(3,1);
  for J = 1:3
    T1(J) = sum(XNODE(1,2:5).*(1+ETA*XNODE(2,2:5)).*ELXY(J,2:5)./4);
    T2(J) = sum(XNODE(2,2:5).*(1+XI *XNODE(1,2:5)).*ELXY(J,2:5)./4);
    XC(J) = sum((1+XI*XNODE(1,2:5)).*(1+ETA*XNODE(2,2:5)).*ELXY(J,2:5)./4);
    XS(J) = ELXY(J,1) - XC(J);
  end
end