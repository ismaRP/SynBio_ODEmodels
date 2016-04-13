function df = odefunc(Ksw1, Ksw2, Kdt11, Kdt12)

%% COMPONENTS
%DNA promoters: Pt1   Pt2    |    Pb1 Pb2
%          RNA: T1    Mt2    |    Mb1 Mb2
%     Proteins:     AT1 Tox2 | Tox1 AT2 Tox3
%   AT Complex:     1AT 1TAT | 2AT 2TAT
%       Growth: g
%     Bacteria: Bac

%% PARAMETERS
%UNITS ARE MINUTES AND MOL. SUPOSSE A VOLUME = 1 UNIT CELL (Arbitrary)
%1 volume cell is supossed to be 2.4*10^-15 liters (average volume of a cell in gro)
%Volume factor to transform M(mol/l) --> Molec/cell
%V=6.022E23 * 2.4E-15;
V=3.612E8;
% INPUT 
%Transcription
aT1 = 1; aMt2 = 1;
%Dissociation constants are given when calling the function
%Translation
bAT1 = 7; bTox2 = 1;
%RNA degradation
T1hl = 6.8 ; Mt2hl = 6.8;
dT1 = log(2)/T1hl ; dMt2 =log(2)/Mt2hl;
%Protein degradation
AT1hl = 30 ; Tox2hl = 60;
dAT1 = log(2)/AT1hl ; dTox2 = log(2)/Tox2hl;

% COMPUTING
%Transcription
aB1 = 1 ; aB2 = 1;
%Dissociation constant of switches are given when calling the function
%Translation
bTox1 = 1 ; bAT2 = 7 ; bTox3 = 1 ;
%RNA degradation
Mb1hl = 6.8 ; Mb2hl = 6.8;
dMb1 = log(2)/Mb1hl ; dMb2 = log(2)/Mb2hl;
%Protein degradation
Tox1hl = 60 ; Tox3hl = 60 ; AT2hl = 30;
dTox1 = log(2)/Tox1hl ; dTox3 = log(2)/Tox3hl ; dAT2 = log(2)/AT2hl;

% TOXIN-ANTITOXIN INTERACTION
kTA1 = 2E6/V ; uTA1 = 7.15E-6; %kTA1 = 2E6/V;
kTAT1 = 2E7/V; uTAT1 = 2.72E-11; %kTAT1 = 2E7/V;
TA1hl = 60; TAT1hl = 60;
dTA1 = log(2)*TA1hl ; dTAT1 = log(2)*TAT1hl;

kTA2 = 2E6/V ; uTA2 = 7.15E-6;
kTAT2 = 2E7/V; uTAT2 = 2.72E-11;
TA2hl = 60; TAT2hl = 60;
dTA2 = log(2)*TA2hl ; dTAT2 = log(2)*TAT2hl;

% BACTERIAL GROWTH AND TOXIN EFFECT
genTime = 40; %min
muMax = log(2)/genTime; dBac = 0.03;
Kt1 = 1E-8*V ; Kt2 = 1E-8*V ; Kt3 = 1E-8*V ; nHill = 2;

df = @(t,y) [
    0; %1 Pt1
    0; %2 Pt2
    0; %3 Pb1
    0; %4 Pb2
    
    aT1  * y(1)      - dT1  * y(5); %5 T1
    aMt2 * y(2)      - dMt2 * y(6); %6 Mt2
    aB1  * (y(5)/(Kdt11+y(5))) * y(3) * Ksw1 + aB1 * (1-Ksw1) * y(3) - dMb1 * y(7); %7 Mb1
    aB2  * (y(5)/(Kdt12+y(5))) * y(4) * Ksw2 + aB2 * (1-Ksw2) * y(4) - dMb2 * y(8); %8 Mb2
    
    %Antitox: +Trad -Degr -TAform +TAbreak
    bAT1 *  y(6)       - dAT1  * y(9)  - kTA1 * y(11) * y(9)  + uTA1 * y(14); %9 AT1
    bAT2 * (y(7)+y(8)) - dAT2  * y(10) - kTA2 * y(12) * y(10) + uTA2 * y(16); %10 AT2
    
    %Tox: +Trad -Degr -TAform +TAbreak -TATform +TATbreak +ATdegrTA +2*ATdegrTAT
    bTox1 * y(7) - dTox1 * y(11) - kTA1 * y(11) * y(9)  + uTA1 * y(14) - kTAT1 * y(11) * y(14) + uTAT1 * y(15) + dAT1 * y(14) + 2 * dAT1 * y(15); %11 Tox1
    bTox2 * y(6) - dTox2 * y(12) - kTA2 * y(12) * y(10) + uTA2 * y(16) - kTAT2 * y(12) * y(16) + uTAT2 * y(17) + dAT2 * y(16) + 2 * dAT2 * y(17); %12 Tox2
    
    bTox3 *  y(8) - dTox3 * y(13)%13 Tox3
    
    %TA:  +TAform  -TAbreak  -ATdegrWithinTA -TATform - TAdegr +TATbreak
    %TAT: +TATform -TATbreak -ATdegrWithinTAT         - TATdegr
    kTA1  * y(11) * y(9)  - uTA1  * y(14) - dAT1 * y(14) - kTAT1 * y(11) * y(14) - dTA1  * y(14) + uTAT1 * y(15); %14 TA1
    kTAT1 * y(14) * y(11) - uTAT1 * y(15) - dAT1 * y(15)                         - dTAT1 * y(15); %15 TAT1
    
    kTA2  * y(12) * y(10) - uTA2  * y(16) - dAT2 * y(16) - kTAT1 * y(12) * y(16) - dTA2  * y(16) + uTAT2 * y(17); %16 TA2
    kTAT2 * y(16) * y(12) - uTAT2 * y(17) - dAT2 * y(17)                         - dTAT2 * y(17); %17 TAT2
    
    %((Kt1^nHill)/((Kt1^nHill) + (y(11)^nHill))) * ((Kt2^nHill)/((Kt2^nHill) + (y(12)^nHill))) * ((Kt3^nHill)/((Kt3^nHill) + (y(13)^nHill))) * muMax * y(18) - dBac * y(18) * (((y(11)+y(12)+y(13))^nHill)/((Kt1^nHill)+(y(11)+y(12)+y(13))^nHill)); %18 Bac
    muMax * y(18) - dBac * y(18) * (((y(11)+y(12)+y(13))^nHill)/((Kt1^nHill)+((y(11)+y(12)+y(13))^nHill))); %18 Bac
    
];

end

