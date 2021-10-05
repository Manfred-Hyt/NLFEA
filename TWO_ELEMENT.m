%
% Two-element example page 131
%
clear
clc
% Nodal coordinates
XYZ=[0 0 0; 1 0 0; 1 1 0; 0 1 0;
    0 0 1; 1 0 1; 1 1 1; 0 1 1;
    0 0 2; 1 0 2; 1 1 2; 0 1 2]*0.01;
%
% Element connectivity
LE=[1 2 3 4 5 6 7 8;
    5 6 7 8 9 10 11 12];
%
% External forces [Node, DOF, Value]
EXTFORCE=[9 3 10.0E3; 10 3 10.0E3; 11 3 10.0E3; 12 3 10.0E3];
%
% Prescribed displacements [Node, DOF, Value]
SDISPT=[1 1 0;1 2 0;1 3 0;2 2 0;2 3 0;3 3 0;4 1 0;4 3 0];
%
% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS=[0.0 0.8 0.4 0.0 0.8; 0.8 1.1 0.1 0.8 1.1]';
%
% Material properties PROP=[LAMDA MU BETA H Y0]
MID=31;
PROP=[110.747E9 80.1938E9 0.0 1.E9 4.0E8];
%
% Set program parameters
ITRA=70; ATOL=1.0E5; NTOL=6; TOL=1E-6;
%
% Calling main function
NOUT = fopen('output.txt','w');
NLFEA(ITRA, TOL, ATOL, NTOL, TIMS, NOUT, MID, PROP, EXTFORCE, SDISPT, XYZ, LE);
fclose(NOUT);