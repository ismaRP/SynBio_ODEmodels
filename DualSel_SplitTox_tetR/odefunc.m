function df = odefunc(kdTF, kdCI, aY)
%% VARIABLES

%DNA promoters P1 P2 P3 Py P4 P5
%RNA M1 M2 M3 My M4 tR
%Prots TF tetR TC TN CI rel B BI Tox
%Bacteria Bac

%% PARAMETERS
%UNITS ARE MINUTES AND MOL. SUPOSSE A VOLUME = 1 UNIT CELL (Arbitrary)
%1 volume cell is supossed to be 2.4*10^-15 liters (average volume of a cell in gro)
%Volume factor to transform M(mol/l) --> Molec/cell
V=6.022E23 * 2.4E-15;

% INPUT PLASMID %

%Transcription
a1 = 5; a2 = 10; a3 = 4;
%Dissociation constant
kdtetR = 0.01;
%Translation
bTF = 0.8; IbTN = 1; IbTC = 1; btetR = 1;
%RNA degradation
M1hl = 6.8; M2hl = 6.8; M3hl = 6.8;
dM1 = log(2)/M1hl; dM2 = log(2)/M2hl; dM3 = log(2)/M3hl;
%Prot degradation
TFhl = 5*60; ToxNhl = 5*60; ToxChl = 5*60; tetRhl = 5*60;
dTF = log(2)/TFhl; dTN = log(2)/ToxNhl; dTC = log(2)/ToxChl; dtetR = log(2)/tetRhl;

% COMPUTING PLASMID %

%Transcription
aYp = 15; a4 = 3; a5 = 3;
% Dissociation constant
kdtR = 0.01;
%Translation
CbTN = 1; CbTC = 1; bCI = 3; brel = 1;
%RNA degradation
Myhl = 6.8; M4hl = 6.8; tRhl = 6.8;
dMy = log(2)/Myhl; dM4 = log(2)/M4hl; dtR = log(2)/tRhl;
%Prot deggradation
relhl = 5*60; CIhl = 5*60;
drel = log(2)/relhl; dCI = log(2)/CIhl;

% Intein reaction
kb = 2.8E7/V; kub = 2.8E5/V;
k1 = 4.7*60E-3; k2 = 7.03*60E-3; k3 = 3.6*60E-4;
Toxhl = 5*60;
dTox = log(2)/Toxhl;

%Bacterial growth and Toxin effect
genTime = 30; %min
muMax = log(2)/genTime;
Kt = 1E-8*V;
nHill = 2;
dBac = 0.0167;

df = @(t,y) [
	0;%1 P1
	0;%2 P2
	0;%3 P3
	0;%4 Py
	0;%5 P4
	0;%6 P5
	
	a1 * y(1)                                                                - dM1 * y(7); %7 M1
	a2 * y(2)                                                                - dM2 * y(8); %8 M2
	a3 * (1-(y(17)/(y(17)+kdCI))) * y(3)                                     - dM3 * y(9); %9 M3
	aYp * (y(13)/(y(13)+kdTF)) * y(4) + aY * (1-(y(13)/(y(13)+kdTF))) * y(4) - dMy * y(10); %10 My
    a4 * (1-(y(14)/(y(14)+kdtetR))) * y(5)                                   - dM4 * y(11); %11 M4
	a5 * (1-(y(17)/(y(17)+kdCI))) * y(6)                                     - dtR * y(12); %12 tR
	
	bTF * y(7)                                                               - dTF * y(13); %13 TF
	btetR * y(8)                                                             - dtetR * y(14); %14 tetR
	IbTC * y(9) + CbTC * y(10) + kub * y(19)                                 - kb * y(15) * y(16) - dTC * y(15); %15 TC
    IbTN * (y(12)/(y(12)+kdtR)) * y(9) + CbTN * y(11) + kub * y(19)          - kb * y(15) * y(16) - dTN * y(16); %16 TN
	bCI * y(10)                                                              - dCI * y(17); %17 CI
	brel * y(10)                                                             - drel * y(18); %18 rel
    
    kb * y(15) * y(16) + k2 * y(20)                                          - kub * y(19) - k1 * y(19); %19 B
    k1 * y(19)                                                               - k2 * y(20) - k3 * y(20); %20 BI
    k3 * y(20)                                                               - dTox * y(21); %21 Tox
    
    ((Kt^nHill)/((Kt^nHill) + (y(21)^nHill))) * muMax * y(22)                - dBac * ((y(21)^nHill)/((Kt^nHill) + (y(21)^nHill))) * y(22); %22 Bac
    % muMax * y(22)                - dBac * ((y(21)^nHill)/((Kt^nHill) + (y(21)^nHill))) * y(22); %22 Bac
    % ((Kt^nHill)/((Kt^nHill) + (y(21)^nHill))) * muMax * y(22) %22 Bac


];
