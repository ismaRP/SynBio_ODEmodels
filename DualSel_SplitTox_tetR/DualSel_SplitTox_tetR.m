%%%ODE SYSTEM FOR DualSel_SplitTox_2 --> tetR modification

%% VARIABLES

%DNA promoters P1 P2 P3 Py P4 P5;
%RNA M1 M2 M3 My M4 tR;
%Prots TF tetR TC TN CI rel B BI Tox

%% PARAMETERS
%UNITS ARE MINUTES AND MOL. SUPOSSE A VOLUME = 1 UNIT CELL (Arbitrary)
%1 volume cell is supossed to be 2.4*10^-15 liters (average volume of a cell in gro)
%Volume factor to transform M(mol/l) --> Molec/cell
V=6.022E23 * 2.4E-15;

%I will add a growth 

iPlasmids = 10;
cPlasmids = 25;
n = 4; %number of different kdTF
r = 2; %number of different kdCI

% INPUT PLASMID %

%Transcription
a1 = 3; a2 = 10; a3 = 2;

%Dissociation constant
if iPlasmids~=0
    kdCI = zeros(r); kdCI = kdCI(1,:);
    kdCI(r)=0.11E-9 * V;
    for i=1:r-1
        kdCI(i)=50*(r-i);
    end
else
    kdCI=0.1;
end
kdtetR = 0.1;
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
if iPlasmids~=0
    kdTF = zeros(n); kdTF = kdTF(1,:);
    kdTF(n) = 0.11E-9*V;
    kdTF(n-1) = 50;
    for i=1:n-2
        kdTF(i) = 100*(n-i-1);
    end
else
    kdTF=0.1;
end
if iPlasmids==0
    alfaY = zeros(n); alfaY = alfaY(1:n);
    alfaY(1)=0.01;
    for i=2:n
        alfaY(i) = alfaY(i-1)/2;
    end
else
    alfaY=0;
end
aY=0;
aYp = 5; a4 = 3; a5 = 3;
% Dissociation constant
kdtR = 0.1;
%Translation
CbTN = 1; CbTC = 1; bCI = 3; brel = 1;
%RNA degradation
Myhl = 6.8; M4hl = 6.8; tRhl = 6.8;
dMy = log(2)/Myhl; dM4 = log(2)/M4hl; dtR = log(2)/tRhl;
%Prot deggradation
relhl = 5*60; CIhl = 5*60;
drel = log(2)/relhl; dCI = log(2)/CIhl;

% Intein reaction
kb = 2.8E7/V;
kub = 2.8E5/V;
k1 = 4.7*60E-3; k2 = 7.03*60E-3; k3 = 3.6*60E-4;
Toxhl = 5*60;
dTox = log(2)/Toxhl;

% P1 = 0; y(1)
% P2 = 0; y(2)
% P3 = 0; y(3)
% Py = 0; y(4)
% P4 = 0; y(5)
% P5 = 0; y(6)

% M1 = a1 * P1 - dM1 * M1;												y(7)
% M2 = a2 * P2 - dM2 * M2;												y(8)
% M3 = a3 * (1-(CI/(CI+kdCI))) * P3 - dM3 * M3;							y(9)
% My = aYp * (TF/(1+kdTF)) * Py + aY * (1-(TF/(1+kdTF))) * Py - dMy * My;y(10)
% M4 = a4 * (1-(tetR/(tetR+kdtetR))) * P4 - dM4 * M4;					y(11)
% tR = a5 * (1-(CI/(CI+kdCI))) * P5 - dtR * tR;							y(12)

% TF = bTF * M1 - dTF * TF;												y(13)
% tetR = btetR * M2 - dtetR * tetR;										y(14)
% TC = IbTC * M3 + CbTC * My - dTC * TC; 								y(15)
% TN = IbTN * (tR/(tR+kdtR)) * M3 + CbTN * M4 - dTN * TN; 				y(16)
% CI = bCI * My - dCI * CI;												y(17)
% rel = brel * My - drel * rel;											y(18)

% B = TC*TN*kb - B * kub                                                y(19)
% BI = B * k1 - BI * k3 - BI * k2                                       y(20)
% Tox = BI * k3 - Tox * dTox                                            y(21)

f = @(t,y) [
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
    IbTN * (y(12)/(y(12)+kdtR)) * y(9) + CbTC * y(11) + kub * y(19)          - kb * y(15) * y(16) - dTN * y(16); %16 TN
	bCI * y(10)                                                              - dCI * y(17); %17 CI
	brel * y(10)                                                             - drel * y(18); %18 rel
    
    kb * y(15) * y(16) + k2 * y(20)                                          - kub * y(19) - k1 * y(19); %19 B
    k1 * y(19)                                                               - k2 * y(20) - k3 * y(20); %20 BI
    k3 * y(20)                                                               - dTox * y(21); %21 Tox
];
y0 = [iPlasmids, iPlasmids, iPlasmids, cPlasmids, cPlasmids, cPlasmids, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

if iPlasmids==0
    disp('--------------INPUT = 0----------------');
    disp(alfaY);
    figure(1);
    xlabel('Time'); ylabel('Tox');title ('Tox vs Time. Input = 0');
    figure(2);
    xlabel('Time'); ylabel('Rel');title ('Rel vs Time. Input = 0');
    hold all;
    for j = 1:length(alfaY)
        disp(alfaY(j));
        aY = alfaY(i);
        [t,y]=ode15s(f, [0 1000], y0);
        figure(1)
        plot (t, y(:,21), 'LineWidth', 1.5); %Toxin
        figure(2)
        plot (t, y(:,18), 'LineWidth', 1.5); %Rel
    end
else
    disp('-------------INPUT = 1----------------');
    figure(1);
    title('Toxin vs Time (All combinations). Input=1');
    xlabel('Time'); ylabel('Toxin');
    [t,y]=ode15s(f, [0 1000], y0);
    plot (t, y(:,21), 'LineWidth', 1.5);   
end




