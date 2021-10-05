%
% One element example
%
% Nodal coordinates
XYZ=[0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
%
% Element connectivity
LE=[1 2 3 4 5 6 7 8];
%
% External forces [Node, DOF, Value]
EXTFORCE=[5 3 10.0E3; 6 3 10.0E3; 7 3 10.0E3; 8 3 10.0E3];
%EXTFORCE=[];
%
% Prescribed displacements [Node, DOF, Value]
SDISPT=[1 1 0;1 2 0;1 3 0;2 2 0;2 3 0;3 3 0;4 1 0;4 3 0];
%SDISPT=[1 1 0;1 2 0;1 3 0;2 2 0;2 3 0;3 3 0;4 1 0;4 3 0;5 3 0.5;6 3 0.5;7 3 0.5;8 3 0.5];
%
% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS=[0.0 0.5 0.1 0.0 0.5; 0.5 1.0 0.1 0.5 1.0]';
%TIMS=[0.0 1.0 1.0 0.0 1.0]';
%
% Contact elements
LEC=[];
%
% Material properties
% MID:0(Linear elastic) 
%        PROP=[LAMBDA NU]
%MID=0;
%PROP=[1.1538E6, 7.6923E5];
%     1(Linear hardening) 2(Saturated hardening)
%     31(Linear hardening finite deformation)
%     32(Saturated hardening finite deformation)
%        PROP=[LAMDA MU BETA H Y0 CE1 FTOL]
MID=1;
PROP=[110.747E3   80.1938E3   0.0    10.E3  50.0E3  1.0  0.0001];
%     -N(Hyperelasticity with N parameters)
%        PROP=[A10 A01 D]
%MID=-2;
%PROP=[80   20   10000];
%
% Set program parameters
ITRA=30; ATOL=1.0E5; NTOL=6; TOL=1E-6;
clc
%
% Calling main function
NOUT = fopen('output.txt','w');
NLFEA(ITRA, TOL, ATOL, NTOL, TIMS, NOUT, MID, PROP, EXTFORCE, SDISPT, XYZ, LE);
fclose(NOUT);