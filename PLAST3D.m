function PLAST3D(MID, PROP, ETAN, UPDATE, LTAN, NE, NDOF, XYZ, LE)
%***********************************************************************
%  MAIN PROGRAM COMPUTING GLOBAL STIFFNESS MATRIX RESIDUAL FORCE FOR
%  PLASTIC MATERIAL MODELS
%***********************************************************************
%%
  global DISPDD DISPTD FORCE GKF SIGMA XQ
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
    DSPD=DISPDD(IDOF);
    DSP=reshape(DSP,NDOF,8);
    DSPD=reshape(DSPD,NDOF,8);
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
      % Previous converged history variables
      if MID > 30
        NALPHA=3; 
        STRESSN=SIGMA(7:12,INTN);
      else
        NALPHA=6; 
        STRESSN=SIGMA(1:6,INTN);
      end
      ALPHAN=XQ(1:NALPHA,INTN);
      EPN=XQ(NALPHA+1,INTN);
      %
      % Strain increment
      if MID == 2 || MID == 31
        F=DSP*SHPD' + eye(3);
        SHPD=inv(F)'*SHPD;
      end
      DEPS=DSPD*SHPD';
      DDEPS=[DEPS(1,1) DEPS(2,2) DEPS(3,3) ...
             DEPS(1,2)+DEPS(2,1) DEPS(2,3)+DEPS(3,2) DEPS(1,3)+DEPS(3,1)]';
      %
      % Computer stress, back stress & effective plastic strain
      if MID == 1
        % Infinitesimal plasticity
        [STRESS, ALPHA, EP]=combHard(PROP,ETAN,DDEPS,STRESSN,ALPHAN,EPN);
      elseif MID == 2
        % Plasticity with finite rotation
        FAC=FAC*det(F);
        [STRESSN, ALPHAN] = rotatedStress(DEPS, STRESSN, ALPHAN);
        [STRESS, ALPHA, EP]=combHard(PROP,ETAN,DDEPS,STRESSN,ALPHAN,EPN);
      elseif MID == 31
        [STRESS, B, ALPHA, EP]=mulPlast(PROP,ETAN,DEPS,STRESSN,ALPHAN,EPN);
      end
      %
      % Update plastic variables
      if UPDATE
        SIGMA(1:6,INTN)=STRESS;
        XQ(:,INTN)= [ALPHA; EP];
        if MID > 30
          SIGMA(7:12,INTN)=B;
        end
        continue;
      end
      %
      % Add residual force and tangent stiffness matrix
      BM=zeros(6,24);
      BG=zeros(9,24);
      for I=1:8
        COL=(I-1)*3+1:(I-1)*3+3;
        BM(:,COL)=[SHPD(1,I) 0         0;
                   0         SHPD(2,I) 0;
                   0         0         SHPD(3,I);
                   SHPD(2,I) SHPD(1,I) 0;
                   0         SHPD(3,I) SHPD(2,I);
                   SHPD(3,I) 0         SHPD(1,I)];
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
      FORCE(IDOF) = FORCE(IDOF) - FAC*BM'*STRESS;
      %
      % Tangent stiffness
      if LTAN
        if MID == 1
          DTAN=combHardTan(PROP,ETAN,DDEPS,STRESSN,ALPHAN,EPN);
          EKF = BM'*DTAN*BM;
        elseif MID == 2
          DTAN=combHardTan(PROP,ETAN,DDEPS,STRESSN,ALPHAN,EPN);
          CTAN=[-STRESS(1) STRESS(1) STRESS(1) -STRESS(4) 0 -STRESS(6);
                STRESS(2) -STRESS(2) STRESS(2) -STRESS(4) -STRESS(5) 0;
                STRESS(3) STRESS(3) -STRESS(3) 0 -STRESS(5) -STRESS(6);
                -STRESS(4) -STRESS(4) 0 -0.5*(STRESS(1)+STRESS(2)) -0.5*STRESS(6) -0.5*STRESS(5);
                0 -STRESS(5) -STRESS(5) -0.5*STRESS(6) -0.5*(STRESS(2)+STRESS(3)) -0.5*STRESS(4);
                -STRESS(6) 0 -STRESS(6) -0.5*STRESS(5) -0.5*STRESS(4) -0.5*(STRESS(1)+STRESS(3))];
          SIG=[STRESS(1) STRESS(4) STRESS(6);
               STRESS(4) STRESS(2) STRESS(5);
               STRESS(6) STRESS(5) STRESS(3)];
          SHEAD=kron(eye(3),SIG);
          EKF = BM'*(DTAN+CTAN)*BM + BG'*SHEAD*BG;
        elseif MID == 31
          DTAN=mulPlastTan(PROP,ETAN,DEPS,STRESSN,ALPHAN,EPN);
          SIG=[STRESS(1) STRESS(4) STRESS(6);
               STRESS(4) STRESS(2) STRESS(5);
               STRESS(6) STRESS(5) STRESS(3)];
          SHEAD=kron(eye(3),SIG);
          %
          EKF = BM'*DTAN*BM + BG'*SHEAD*BG;
        end
        GKF(IDOF,IDOF)=GKF(IDOF,IDOF)+FAC*EKF;
      end
    end,        end,        end
  end
end