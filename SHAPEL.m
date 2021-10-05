function [SF, GDSF, DET] = SHAPEL(XI, ELXY)
%*************************************************************************
% Compute shape function, derivatives, and determinant of hexahedron element
%*************************************************************************
%%
  XNODE=[-1  1  1 -1 -1  1  1 -1;
         -1 -1  1  1 -1 -1  1  1;
         -1 -1 -1 -1  1  1  1  1];
  QUAR = 0.125;
  SF=zeros(8,1);
  DSF=zeros(3,8);
  for I=1:8
    XP = XNODE(1,I);
    YP = XNODE(2,I);
    ZP = XNODE(3,I);
    %
    XI0 = [1+XI(1)*XP 1+XI(2)*YP 1+XI(3)*ZP];
    %
    SF(I) = QUAR*XI0(1)*XI0(2)*XI0(3);
    DSF(1,I) = QUAR*XP*XI0(2)*XI0(3);
    DSF(2,I) = QUAR*YP*XI0(1)*XI0(3);
    DSF(3,I) = QUAR*ZP*XI0(1)*XI0(2);
  end
  GJ = DSF*ELXY;
  DET = det(GJ);
  GJINV=inv(GJ);
  GDSF=GJINV*DSF;
  
end