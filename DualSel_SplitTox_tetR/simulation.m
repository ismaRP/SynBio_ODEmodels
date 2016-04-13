%% Simulation & Plots
P1=1; P2=2; P3=3; Py=4; P4=5; P5=6;
M1=7; M2=8; M3=9; My=10; M4=11; tr=12;
TF=13; tetR=14; TC=15; TN=16; CI=17; rel=18; B=19; BI=20; Tox=21;
Bac=22;

r=8;
n=8;
iPlasmids=10;
cPlasmids=25;
initialBac = 1000;
V=6.022E23 * 2.4E-15;

if iPlasmids~=0
    kdCI = zeros(r); kdCI = kdCI(1,:);
    kdCI(r)=0.11E-9 * V;
    for i=1:r-1
        kdCI(i)=50*(r-i);
    end
else
    kdCI=0.1;
end
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
    aY = zeros(n); aY = aY(1:n);
    aY(1)=0.1;
    for i=2:n
        aY(i) = aY(i-1)/2;
    end
    % aY=[0.01 0.001];
else
    aY=0;
end

%Set initial Bac depending on the gate
y0 = [iPlasmids, iPlasmids, iPlasmids, cPlasmids, cPlasmids, cPlasmids, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, initialBac];

%%SIMULATIONS
tspan = [0 500];
if iPlasmids==0
    disp('--------------INPUT = 0----------------');
    disp(aY);
    figure(1);
    figure(2);
    figure(3);
    for j = 1:length(aY)
        disp(aY(j));
        deqs=odefunc(kdTF, kdCI, aY(j));
        [t,y]=ode15s(deqs, tspan, y0);
        figure(1)
        plot (t, y(:,Tox), 'LineWidth', 1.5);hold on; %Toxin
        figure(2)
        plot (t, y(:,rel), 'LineWidth', 1.5);hold on; %Rel
        figure(3)
        plot (t,y(:,Bac), 'LineWidth', 1.5);hold on; %Bac
    end
    figure(1);
    xlabel('Time'); ylabel('Tox'); title ('Tox vs Time. Input = 0');
    figure(2);
    xlabel('Time'); ylabel('Rel'); title ('Rel vs Time. Input = 0');
    figure(3);
    xlabel('Time'); ylabel('Bac'); title ('Bac vs Time. Input = 0');
else
    disp('-------------INPUT = 1----------------');
    figure(1);
    figure(2);
    figure(3);
    figure(4);
    figure(5);
    figure(6);
    figure(9);
    MAXtoxin = zeros(r,n);
    for j = 1:length(kdCI)
        for i = 1:length(kdTF)
            disp(strcat('kdCI= ', num2str(kdCI(j))));
            disp(strcat('kdTF= ', num2str(kdTF(i))));
            deqs=odefunc(kdTF(i), kdCI(j), aY);
            [t,y]=ode15s(deqs, tspan, y0);
            
            figure(1)
            plot(t,y(:,Tox), 'LineWidth', 1.5); hold on;
            
            figure(9)
            plot(t,y(:,Bac), 'LineWidth', 1.5); hold on;
            
            MAXtoxin(j,i)=max(y(:,Tox));
            disp(MAXtoxin);
            if j==r %Plots for the best kdCI
                figure(2)
                plot(t, y(:,Tox), 'LineWidth', 1.5);hold on;
                
                TFval = y(:,TF);
                PyP = cPlasmids*(TFval ./ (TFval+kdTF(i)));
                figure(3)
                plot(t, PyP, 'LineWidth', 1.5);hold on;
                figure(4)
                plot(TFval, PyP, 'LineWidth', 1.5);hold on;
            end
        end %Plots for the best kdTF
        figure(5)
        plot(y(:,TF), y(:,Tox), 'LineWidth', 1.5);hold on;
        figure(6)
        plot(t, y(:,Tox), 'LineWidth', 1.5); hold on;
    end
    
    figure(1);
    title('Toxin vs Time (All combinations). Input=1');
    xlabel('Time'); ylabel('Toxin');
    figure(2);
    title('Toxin vs Time (Best CI). Input=1');
    xlabel('Time'); ylabel('Toxin');
    figure(3);
    title('PyP vs Time (Best CI). Input=1')
    xlabel('Time'); ylabel('PyP');
    figure(4)
    title('PyP vs TF (Best CI). Input=1')
    xlabel('TF'); ylabel('PyP');
    figure(5)
    title('Toxin vs TF (Best TF). Input=1');
    xlabel('TF'); ylabel('Toxin');
    figure(6)
    title('Toxin vs Time (Best TF). Input=1');
    xlabel('Time'); ylabel('Toxin');
    figure(7)
    surf(kdTF,kdCI,MAXtoxin);
    title('Toxin peak VS kdCI and kdTF');
    xlabel('kdTF'); ylabel('kdCI');
    figure(8)
    plot(t,y(:,TF), 'LineWidth', 1.5);
    title('TF vs Time');
    xlabel('Time'); ylabel('TF');
    figure(9)
    xlabel('Time'); ylabel('Bac'); title ('Bac vs Time. Input = 0');
end
