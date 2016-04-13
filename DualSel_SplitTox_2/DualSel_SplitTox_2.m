%% DUAL SELECTION SPLIT TOX CIRCUIT

%This version models translation with a general translation rate
%This version models TF binding with a Kd, instead of sigma and sigma
%reverse

clear;
%% this block fetches the framework: do not edit
urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
    'part_composition_setup.m');
rehash;
part_composition_setup('v5');

%% DECLARATION 
sys = BioSystem();

%% VARIABLES
%UNITS ARE MINUTES AND MOL. SUPOSSE A VOLUME = 1 UNIT CELL (Arbitrary)
%1 volume cell is supossed to be 2.4*10^-15 liters (average volume of a cell in gro)
%Volume factor to transform M(mol/l) --> Molec/cell
V=6.022E23 * 2.4E-15;

%I will add a growth 

iPlasmids = 0;
cPlasmids = 25;
n = 2; %number of different kdTF
r = 2; %number of different kdCI

% Input plasmid %

M1Length = 1000; tR1Length = 50; M31Length = 1000; M32Length = 2000; mediaI = (M31Length+M32Length)/2;
TFLength = 300; ToxCLength = 300; ToxNLength = 300;
%Transcription
alfa1 = 3*60/M1Length; alfa2 = 4*60/tR1Length; alfa3 = 15*60/mediaI;

%Dissociation constant
kdCI = zeros(r); kdCI = kdCI(1,:);
kdCI(r)=0.11E-9 * V;
for i=1:r-1
    kdCI(i)=50*(r-i);
end

kdtR2 = 1;
%Translation 
betaTF = 0.8; IbetaToxN = 1; IbetaToxC = 1;
%RNA degradation
M1hl = 6.8; M31hl = 6.8; M32hl = 6.8; tR1hl = 6.8;
deltaM1 = log(2)/M1hl; deltaM31 = log(2)/M31hl; deltaM32 = log(2)/M32hl; deltatR1 = log(2)/tR1hl;
%Prot degradation
TFhl = 5*60; ToxNhl = 5*60; ToxChl = 5*60;
deltaTF = log(2)/TFhl; deltaToxN = log(2)/ToxNhl; deltaToxC = log(2)/ToxChl;

% Computing plasmid %

My1Length = 3000; My2Length = 4000; mediaC = (My1Length+My2Length)/2;
tR2Length = 50;
CILength = 300; RelLength = 300;
%Transcription
kdTF = zeros(n); kdTF = kdTF(1,:);
kdTF(n) = 0.11E-9*V;
kdTF(n-1) = 50;
for i=1:n-2
    kdTF(i) = 100*(n-i-1);
end
alfaY = zeros(n); alfaY = alfaY(1:n);
alfaY(1)=0.01;
for i=2:n
    alfaY(i) = alfaY(i-1)/2;
end

alfaYp = 20*60/mediaC; alfa4 = 4*60/tR2Length;
% Dissociation constant
kdtR1 = 0.1;
%Translation
CbetaToxN = 1; CbetaToxC = 1; betaCI = 3; betaRel = 1;
%RNA degradation
My1hl = 6.8; My2hl = 6.8; tR2hl = 6.8;
deltaMy1 = log(2)/My1hl; deltaMy2 = log(2)/My2hl; deltatR2 = log(2)/tR2hl;
%Prot deggradation
Relhl = 5*60; CIhl = 5*60;
deltaRel = log(2)/Relhl; deltaCI = log(2)/CIhl;

% Intein reaction
kb = 2.8E7/V;
k1 = 4.7*60E-3; k2 = 7.03*60E-3; k3 = 3.6*60E-4;
Toxhl = 5*60;
deltaTox= log(2)/Toxhl;

sys.AddConstant('alfa1', alfa1), sys.AddConstant('alfa2', alfa2), sys.AddConstant('alfa3', alfa3);
sys.AddConstant('kdCI', kdCI(1)), sys.AddConstant('kdtR2', kdtR2);
sys.AddConstant('betaTF', betaTF), sys.AddConstant('IbetaToxN', IbetaToxN), sys.AddConstant('IbetaToxC', IbetaToxC);
sys.AddConstant('deltaM1', deltaM1), sys.AddConstant('deltaM31', deltaM31), sys.AddConstant('deltaM32', deltaM32), sys.AddConstant('deltatR1', deltatR1);
sys.AddConstant('deltaTF', deltaTF), sys.AddConstant('deltaToxN', deltaToxN), sys.AddConstant('deltaToxC', deltaToxC);

sys.AddConstant('alfaY', alfaY(1)), sys.AddConstant('alfaYp', alfaYp);
sys.AddConstant('alfa4', alfa4);
sys.AddConstant('kdtR1', kdtR1); sys.AddConstant('kdTF', kdTF(1));
sys.AddConstant('CbetaToxN', CbetaToxN), sys.AddConstant('CbetaToxC', CbetaToxC), sys.AddConstant('betaCI', betaCI), sys.AddConstant('betaRel', betaRel);
sys.AddConstant('deltaMy1', deltaMy1), sys.AddConstant('deltaMy2', deltaMy2), sys.AddConstant('deltatR2', deltatR2);
sys.AddConstant('deltaRel', deltaRel), sys.AddConstant('deltaCI', deltaCI);

sys.AddConstant('kb', kb);
sys.AddConstant('k1', k1), sys.AddConstant('k2', k2), sys.AddConstant('k3', k3);
sys.AddConstant('deltaTox', deltaTox);
%% COMPOSITORS
dP1dt = sys.AddCompositor('P1', iPlasmids); dP2dt = sys.AddCompositor('P2', iPlasmids); dP3dt = sys.AddCompositor('P3', iPlasmids);

dM1dt = sys.AddCompositor('M1', 0);
dtR1dt = sys.AddCompositor('tR1', 0);

dM31dt = sys.AddCompositor('M31', 0);
dM32dt = sys.AddCompositor('M32', 0);

dTFdt = sys.AddCompositor('TF', 0);
dToxNdt = sys.AddCompositor('ToxN', 0);
dToxCdt = sys.AddCompositor('ToxC', 0);

dPydt = sys.AddCompositor('Py', cPlasmids);
dP4dt = sys.AddCompositor('P4', cPlasmids);

dMy1dt = sys.AddCompositor('My1', 0);

dMy2dt = sys.AddCompositor('My2', 0);

dtR2dt = sys.AddCompositor('tR2', 0);
dReldt = sys.AddCompositor('Rel', 0);
dCIdt = sys.AddCompositor('CI', 0);

dBdt = sys.AddCompositor('B', 0);
dBIdt = sys.AddCompositor('BI', 0);
dToxdt = sys.AddCompositor('Tox', 0);

dDegrdt = sys.AddCompositor('Degr', 0);

%% PARTS
%Input Transcrition
P1Transcr = Part('P1 -alfa1> M1 + P1', [dP1dt dM1dt], ... 
    [Rate('0') Rate('alfa1 * P1')]);
sys.AddPart(P1Transcr);

P2Transcr = Part('P2 -alfa2> tR1 + P2', [dP2dt dtR1dt], ...
    [Rate('0') Rate('alfa2 * P2')]);
sys.AddPart(P2Transcr);

P3Transcr1 = Part('P3 -alfa3> M31 + P3', [dP3dt dM31dt], ...
    [Rate('0') Rate('alfa3 * (1-(CI/(kdCI+CI))) * (1-(tR2/(kdtR2+tR2))) * P3')]);
sys.AddPart(P3Transcr1);

P3Transcr2 = Part('P3 -alfa3> M32 + P3', [dP3dt dM32dt], ...
    [Rate('0') Rate('alfa3 * (1-(CI/(kdCI+CI))) * (tR2/(kdtR2+tR2)) * P3')]);
sys.AddPart(P3Transcr2);

%Input Translation
M1elong = Part('M1 -betaTF> M1 + TF', [dM1dt dTFdt], ...
    [Rate('0') Rate('betaTF * M1')]);
sys.AddPart(M1elong);

M31elong = Part('M31 -IbetaToxC> M31 + ToxC', [dM31dt dToxCdt], ...
    [Rate('0') Rate('IbetaToxC * M31')]);
sys.AddPart(M31elong);

M32Celong = Part('M32 -IbetaToxC> M32 + ToxC', [dM32dt dToxCdt], ...
    [Rate('0') Rate('IbetaToxC * M32')]);
sys.AddPart(M32Celong);

M32Nelong = Part('M32 -IbetaToxN> M32 + ToxN', [dM32dt dToxNdt], ...
    [Rate('0') Rate('IbetaToxN * M32')]);
sys.AddPart(M32Nelong);

%Input RNA degr
M1Degr = Part('M1 -deltaM1> Degr', [dM1dt dDegrdt], ...
    [Rate('-deltaM1 * M1') Rate('0')]);
sys.AddPart(M1Degr);
tR1Degr = Part('tR1 -deltatR1> Degr', [dtR1dt dDegrdt], ...
    [Rate('-deltatR1 * tR1') Rate('0')]);
sys.AddPart(tR1Degr);
M31Degr = Part('M31 -deltaM31> Degr', [dM31dt dDegrdt], ...
    [Rate('-deltaM31 * M31') Rate('0')]);
sys.AddPart(M31Degr);
M32Degr = Part('M32 -deltaM32> Degr', [dM32dt dDegrdt], ...
    [Rate('-deltaM32 * M32') Rate('0')]);
sys.AddPart(M32Degr);

%Input Prot degr
TFDegr = Part('TF -deltaTF> Degr', [dTFdt dDegrdt], ...
    [Rate('-deltaTF * TF') Rate('0')]);
sys.AddPart(TFDegr);
ToxCDegr = Part('ToxC -deltaToxC> Degr', [dToxCdt dDegrdt], ...
    [Rate('-deltaToxC * ToxC') Rate('0')]);
sys.AddPart(ToxCDegr);
ToxNDegr = Part('ToxN -deltaToxN> Degr', [dToxNdt dDegrdt], ...
    [Rate('-deltaToxN * ToxN') Rate('0')]);
sys.AddPart(ToxNDegr);

%Py Transcription
if iPlasmids == 0
    PyTranscr1 = Part('Py -alfaY> My1 + Py', [dPydt dMy1dt], ...
        [Rate('0') Rate('Py * alfaY * (tR1/(kdtR1+tR1))')]);
    sys.AddPart(PyTranscr1);
    
    PyTranscr2 = Part('Py -alfaY> My2 + Py', [dPydt dMy2dt], ...
        [Rate('0') Rate('Py * alfaY * (1-(tR1/(kdtR1+tR1)))')]);
    sys.AddPart(PyTranscr2);
else
    PyTranscr1 = Part('Py -alfaYp> My1 + Py', [dPydt dMy1dt], ...
        [Rate('0') Rate('Py * alfaYp * (TF/(kdTF+TF)) * (tR1/(kdtR1+tR1))')]);
    sys.AddPart(PyTranscr1);

    PyTranscr2 = Part('Py -alfaYp> My2 + Py', [dPydt dMy2dt], ...
        [Rate('0') Rate('Py * alfaYp * (TF/(kdTF+TF)) * (1-(tR1/(kdtR1+tR1)))')]);
    sys.AddPart(PyTranscr2);
end

%P4 Transcription
P4Transcr = Part('P4 -alfa4> tR2 + P4', [dP4dt dtR2dt], ...
    [Rate('0') Rate('P4 * alfa4 * (1-(CI/(kdCI + CI)))')]);
sys.AddPart(P4Transcr);

%Translation
%1
My1Celong = Part('My1 -CbetaToxC> My1 + ToxC', [dMy1dt dToxCdt], ...
    [Rate('0') Rate('My1 * CbetaToxC')]);
sys.AddPart(My1Celong);

My1CIelong = Part('My1 -betaCI> My1 + CI', [dMy1dt dCIdt], ...
    [Rate('0') Rate('My1 * betaCI')]);
sys.AddPart(My1CIelong);

My1Relong = Part('My1 -betaRel> My1 + Rel', [dMy1dt dReldt], ...
    [Rate('0') Rate('My1 * betaRel')]);
sys.AddPart(My1Relong);
%2
My2Celong = Part('My2 -CbetaToxC> My2 + ToxC', [dMy2dt dToxCdt], ...
    [Rate('0') Rate('My2 * CbetaToxC')]);
sys.AddPart(My2Celong);

My2CIelong = Part('My2 -betaCI> My2 + CI', [dMy2dt dCIdt], ...
    [Rate('0') Rate('My2 * betaCI')]);
sys.AddPart(My2CIelong);

My2Relong = Part('My2 -betaRel> My2 + Rel', [dMy2dt dReldt], ...
    [Rate('0') Rate('My2 * betaRel')]);
sys.AddPart(My2Relong);

My2Nelong = Part('My2 -CbetaToxN> My2 + ToxN', [dMy2dt dToxNdt], ...
    [Rate('0') Rate('My2 * CbetaToxN')]);
sys.AddPart(My2Nelong);

%RNA Degradation
My1Degr = Part('My1 -deltaMy1> Degr', [dMy1dt dDegrdt], ...
    [Rate('-deltaMy1 * My1') Rate('0')]);
sys.AddPart(My1Degr);

My2Degr = Part('My2 -deltaMy2> Degr', [dMy2dt dDegrdt], ...
    [Rate('-deltaMy2 * My2') Rate('0')]);
sys.AddPart(My2Degr);

tR2Degr = Part('tR2 -deltatR2> Degr', [dtR2dt dDegrdt], ...
    [Rate('-deltatR2 * tR2') Rate('0')]);
sys.AddPart(tR2Degr);

%Prot Degradation
%ToxN and ToxC done
CIDegr = Part('CI -deltaCI> Degr', [dCIdt dDegrdt], ...
    [Rate('-deltaCI * CI') Rate('0')]);
sys.AddPart(CIDegr);

RelDegr = Part('Rel -deltaRel> Degr', [dReldt dDegrdt], ...
    [Rate('-deltaRel * Rel') Rate('0')]);
sys.AddPart(RelDegr);

%Intein reaction
IntBinding = Part('ToxC + ToxN -kb> B', [dToxNdt dToxCdt dBdt], ...
    [Rate('-kb * ToxC * ToxN') Rate('-kb * ToxC * ToxN') Rate('kb * ToxC * ToxN')]);
sys.AddPart(IntBinding);

BIReaction = Part('B -k1> BI', [dBdt dBIdt], ...
    [Rate('-k1 * B') Rate('k1 * B')]);
sys.AddPart(BIReaction);

BIReverse = Part('BI -k2> B', [dBIdt dBdt], ...
    [Rate('-k2 * BI') Rate('k2 * BI')]);
sys.AddPart(BIReverse);

BIResolution = Part('BI -k3> Tox', [dBIdt dToxdt], ...
    [Rate('-k3 * BI') Rate('k3 * BI')]);
sys.AddPart(BIResolution);

ToxDegr = Part('Tox -deltaTox> Degr', [dToxdt dDegrdt], ...
    [Rate('-deltaTox * Tox') Rate('0')]);
sys.AddPart(ToxDegr);


%% RUNNING

tspan=[0 2000];
if iPlasmids == 0
    disp('----------INPUT = 0------------');
    disp(alfaY);
    figure(1);
    figure(2);
    figure(3);
    hold all;
    for i = 1:length(alfaY)
        disp(i);
        sys.ChangeConstantValue('alfaY', alfaY(i));
        [T, Y] = sys.run(tspan);
        figure(1)
        plot (T, Y(:,sys.CompositorIndex('Tox')), 'LineWidth', 1.5);
        hold on;
        figure(2)
        plot (T, Y(:,sys.CompositorIndex('Rel')), 'LineWidth', 1.5);
        hold on;
        figure(3)
        plot (T, Y(:,sys.CompositorIndex('My2')), 'LineWidth', 1.5);
        hold on;
    end
%     legend (strcat('\alpha_y=',num2str(alfaY(1))), strcat('\alpha_y=',num2str(alfaY(2))), strcat('\alpha_y=',num2str(alfaY(3))), strcat('\alpha_y=',num2str(alfaY(4))), strcat('\alpha_y=',num2str(alfaY(5))), strcat('\alpha_y=',num2str(alfaY(6))), strcat('\alpha_y=',num2str(alfaY(7))));
    figure(1)
    xlabel('Time'); ylabel('Tox');title ('Tox vs Time. Input = 0');
    
    figure(2)
    xlabel('Time'); ylabel('Rel');title ('Rel vs Time. Input = 0');
    
    figure(3)
    xlabel('Time'); ylabel('My2'); title ('My2 vs Time. Input = 0');
else
    disp('----------INPUT = 1------------')
    figure(1);
    figure(2);
    figure(3);
    figure(4);
    figure(5);
    figure(6);
    MAXtoxin = zeros(r,n);
    for j = 1:length(kdCI)
        sys.ChangeConstantValue('kdCI', kdCI(j));
        for i = 1:length(kdTF)
            disp(strcat('kdCI= ', num2str(kdCI(j))));
            disp(strcat('kdTF= ', num2str(kdTF(i))));
            sys.ChangeConstantValue('kdTF', kdTF(i));
            [T, Y] = sys.run(tspan);
            figure(1)
            plot (T, Y(:,sys.CompositorIndex('Tox')), 'LineWidth', 1.5);
            hold on;
            MAXtoxin(j,i) = max(Y(:,sys.CompositorIndex('Tox')));
            disp(MAXtoxin);
            if j==r %Plots for the best kdCI
                figure(5)
                plot(T, Y(:,sys.CompositorIndex('Tox')), 'LineWidth', 1.5);
                hold on;
                
                TFval = Y(:,sys.CompositorIndex('TF'));
                PyP = cPlasmids*(TFval ./ (TFval+kdTF(i)));
                figure(6)
                plot(T, PyP, 'LineWidth', 1.5);
                hold on;
                figure(7)
                plot(TFval, PyP, 'LineWidth', 1.5);
                hold on;
            end
        end
        %Plots for the best kdTF
        figure(3)
        plot (Y(:,sys.CompositorIndex('TF')), Y(:,sys.CompositorIndex('Tox')), 'LineWidth', 1.5);
        hold on;
        figure(4)
        plot (T, Y(:,sys.CompositorIndex('Tox')), 'LineWidth', 1.5);
        hold on;
    end
    hold off;
%     X = zeros(n);
%     Y = zeros(n);
%     for i=1:n;
%         X(i,:) = kdTF;
%         Y(:,i) = kdCI;
%     end
    figure(1)
    title('Toxin vs Time (All combinations). Input=1');
    xlabel('Time'); ylabel('Toxin');
    
    figure(2)
    surf(kdTF,kdCI,MAXtoxin);
    title('Toxin peak VS kdCI and kdTF');
    xlabel('kdTF'); ylabel('kdCI');
    
    figure(3)
    title('Toxin vs TF (Best TF). Input=1');
    xlabel('TF'); ylabel('Toxin');
    
    figure(4)
    title('Toxin vs Time (Best TF). Input=1');
    xlabel('Time'); ylabel('Toxin');
    
    figure(5)
    title('Toxin vs Time (Best CI). Input=1');
    xlabel('Time'); ylabel('Toxin');
    
    figure(6)
    title('PyP vs Time (Best CI). Input=1')
    xlabel('Time'); ylabel('PyP');
    
    figure(7)
    title('PyP vs TF (Best CI). Input=1')
    xlabel('TF'); ylabel('PyP');
    
    figure(8)
    plot(T,Y(:,sys.CompositorIndex('TF')), 'LineWidth', 1.5);
    title('TF vs Time');
    xlabel('Time'); ylabel('TF');
    
end
