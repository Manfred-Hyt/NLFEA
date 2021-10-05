%
% Bending of five elements: linear elastic
%
% Nodal coordinates (mm)
XYZ=[0 0 0 ;1 0 0 ;1 1 0 ;0 1 0 ;
     0 0 2 ;1 0 2 ;1 1 2 ;0 1 2 ;
     0 0 4 ;1 0 4 ;1 1 4 ;0 1 4 ;
     0 0 6 ;1 0 6 ;1 1 6 ;0 1 6 ;
     0 0 7 ;1 0 7 ;1 1 8 ;0 1 8 ;
     0 0 10;1 0 10;1 1 10;0 1 10];
%
% Element connectivity
LE=[ 1  2  3  4  5  6  7  8;
     5  6  7  8  9 10 11 12;
     9 10 11 12 13 14 15 16;
    13 14 15 16 17 18 19 20;
    17 18 19 20 21 22 23 24];
%
% External forces [Node, DOF, Value]
EXTFORCE=[21 2 1.0E5; 22 2 1.0E5; 23 2 1.0E5; 24 2 1.0E5];
%
% Prescribed displacements [Node, DOF, Value]
SDISPT=[1 1 0;1 2 0;1 3 0;
        2 1 0;2 2 0;2 3 0;
        3 1 0;3 2 0;3 3 0;
        4 1 0;4 2 0;4 3 0];
%
% Material properties
% MID:0(Linear elastic) PROP=[LAMBDA MNU]
%MID=0;
MID=31;
%
%for elastic(MID=0)
%E=2E11; NU=0.3; LAMBDA=E*NU/((1+NU)*(1-2*NU)); MU=E/(2*(1+NU));
%PROP=[LAMBDA MU];
%for elasto plastic mutiply(MID=31)
Young = 24000; nu=0.2; mu=Young/2/(1+nu); lambda=nu*Young/((1+nu)*(1-2*nu));
beta = 0; H = 1000; sY = 200*sqrt(3);
mp = [lambda mu beta H sY];
%
% Load increments [Start End Increment InitialFactor FinalFactor]
%TIMS=[0.0 1.0 1.0 0.0 1.0]';
TIMS=[0.0 1.0 0.1 0.1 0.2]';
%
% Set program parameters
ITRA=30; ATOL=1.0E5; NTOL=6; TOL=1E-6;
%
% Calling main function
NOUT = fopen('output.txt','w');
NLFEA(ITRA, TOL, ATOL, NTOL, TIMS, NOUT, MID, PROP, EXTFORCE, SDISPT, XYZ, LE);
fclose(NOUT);