function HYPER3D(MID, PROP, UPDATE, LTAN, NE, NDOF, XYZ, LE)
%***********************************************************************
%  MAIN PROGRAM COMPUTING GLOBAL STIFFNESS MATRIX RESIDUAL FORCE FOR
%  PLASTIC MATERIAL MODELS
%***********************************************************************
%%
  global DISPTD FORCE GKF SIGMA
  %
  % Integration points and weights
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
      % Deformation gradient
      F=DSP*SHPD' + eye(3);
      %
      % Computer stress and tangent stiffness
      [STRESS DTAN] = Mooney(F, PROP(1), PROP(2), PROP(3), LTAN);
      %
      % Update plastic variables
      if UPDATE
        SIGMA(:,INTN)=STRESS;
        continue;
      end
      %
      % Add residual force and tangent stiffness matrix
      BN=zeros(6,24);
      BG=zeros(9,24);
      for I=1:8
        COL=(I-1)*3+1:(I-1)*3+3;
        BN(:,COL)=[SHPD(1,I)*F(1,1) SHPD(1,I)*F(2,1) SHPD(1,I)*F(3,1);
                   SHPD(2,I)*F(1,2) SHPD(2,I)*F(2,2) SHPD(2,I)*F(3,2);
                   SHPD(3,I)*F(1,3) SHPD(3,I)*F(2,3) SHPD(3,I)*F(3,3);
                   SHPD(1,I)*F(1,2)+SHPD(2,I)*F(1,1) SHPD(1,I)*F(2,2)+SHPD(2,I)*F(2,1) SHPD(1,I)*F(3,2)+SHPD(2,I)*F(3,1);
                   SHPD(2,I)*F(1,3)+SHPD(3,I)*F(1,2) SHPD(2,I)*F(2,3)+SHPD(3,I)*F(2,2) SHPD(2,I)*F(3,3)+SHPD(3,I)*F(3,2);
                   SHPD(1,I)*F(1,3)+SHPD(3,I)*F(1,1) SHPD(1,I)*F(2,3)+SHPD(3,I)*F(2,1) SHPD(1,I)*F(3,3)+SHPD(3,I)*F(3,1)];
        %
        BG(:,COL)=[SHPD(1,I) 0         0;
                   SHPD(2,I) 0         0;
                   SHPD(3,I) 0         0;
                   0         SHPD(1,I) 0;
                   0         SHPD(2,I) 0;
                   0         SHPD(3,I) 0;
                   0         0         SHPD(1,I);
                   0         0         SHPD(2,I);
                   0         0         SHPD(3,I)];
      end
      %
      % Residual forces
      FORCE(IDOF) = FORCE(IDOF) - FAC*BN'*STRESS;
      %
      % Tangent stiffness
      if LTAN
        SIG=[STRESS(1) STRESS(4) STRESS(6);
             STRESS(4) STRESS(2) STRESS(5);
             STRESS(6) STRESS(5) STRESS(3)];
        SHEAD=zeros(9);
        SHEAD(1:3,1:3)=SIG;
        SHEAD(4:6,4:6)=SIG;
        SHEAD(7:9,7:9)=SIG;
        %
        EKF = BN'*DTAN*BN + BG'*SHEAD*BG;
        GKF(IDOF,IDOF)=GKF(IDOF,IDOF)+FAC*EKF;
      end
    end; end; end;
  end
end