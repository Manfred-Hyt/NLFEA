function ELAST3D(ETAN, UPDATE, LTAN, NE, NDOF, XYZ, LE)
%***********************************************************************
%  MAIN PROGRAM COMPUTING GLOBAL STIFFNESS MATRIX RESIDUAL FORCE FOR
%  PLASTIC MATERIAL MODELS
%***********************************************************************
%%
  global DISPTD FORCE GKF SIGMA
  %
  % Integration points and weights (2-point integration)
  XG=[-0.57735026918963D0, 0.57735026918963D0];
  WGT=[1.00000000000000D0, 1.00000000000000D0];
  %
  % Index for history variables (each integration pt)
  INTN=0;
  %
  %LOOP OVER ELEMENTS, THIS IS MAIN LOOP TO COMPUTE K AND F
  for IE=1:NE
    % Nodal coordinates and incremental displacements
    ELXY=XYZ(LE(IE,:),:);
    % Local to global mapping
    IDOF=zeros(1,24);
    for I=1:8
      II=(I-1)*NDOF+1;
      IDOF(II:II+2)=(LE(IE,I)-1)*NDOF+1:(LE(IE,I)-1)*NDOF+3;
    end
    DSP=DISPTD(IDOF);
    DSP=reshape(DSP,NDOF,8);
    %
    %LOOP OVER INTEGRATION POINTS
    for LX=1:2, for LY=1:2, for LZ=1:2
      E1=XG(LX); E2=XG(LY); E3=XG(LZ);
      INTN = INTN + 1;
      %
      % Determinant and shape function derivatives
      [~, SHPD, DET] = SHAPEL([E1 E2 E3], ELXY);
      FAC=WGT(LX)*WGT(LY)*WGT(LZ)*DET;
      %
      % Strain
      DEPS=DSP*SHPD';
      DDEPS=[DEPS(1,1) DEPS(2,2) DEPS(3,3) ...
             DEPS(1,2)+DEPS(2,1) DEPS(2,3)+DEPS(3,2) DEPS(1,3)+DEPS(3,1)]';
      %
      % Stress
      STRESS = ETAN*DDEPS;
      %
      % Update stress
      if UPDATE
        SIGMA(:,INTN)=STRESS;
        continue;
      end
      %
      % Add residual force and tangent stiffness matrix
      BM=zeros(6,24);
      for I=1:8
        COL=(I-1)*3+1:(I-1)*3+3;
        BM(:,COL)=[SHPD(1,I) 0         0;
                   0         SHPD(2,I) 0;
                   0         0         SHPD(3,I);
                   SHPD(2,I) SHPD(1,I) 0;
                   0         SHPD(3,I) SHPD(2,I);
                   SHPD(3,I) 0         SHPD(1,I)];
      end
      %
      % Residual forces
      FORCE(IDOF) = FORCE(IDOF) - FAC*BM'*STRESS;
      %
      % Tangent stiffness
      if LTAN
        EKF = BM'*ETAN*BM;
        GKF(IDOF,IDOF)=GKF(IDOF,IDOF)+FAC*EKF;
      end
    end, end, end
  end
end