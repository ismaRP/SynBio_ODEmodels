%% COMPONENTS
Pt1=1; Pt2=2; Pb1=3; Pb2=4;
T1=5; Mt2=6; Mb1=7; Mb2=8;
AT1=9; AT2=10; Tox1=11; Tox2=12; Tox3=13;
TA1=14; TAT1=15; TA2=16; TAT2=17;
Bac=18;

n=2; %r=;
iPlasmids=0; cPlasmids=10;
initialBac=1000;
%V=6.022E23 * 2.4E-15;
V=3.612E8;

if iPlasmids~=0 %input=1
    Kdt11=zeros(1,n);
    Kdt12=zeros(1,n);
    Ksw1=1; Ksw2=1;
    Kdt11(1)=0.11E-9*V; Kdt11(n)=1000000;
    Kdt12(1)=0.11E-9*V; Kdt12(n)=1000000;
else %input=0
    Ksw1=zeros(1,n);
    Ksw2=zeros(1,n);
    Kdt11=1; Kdt12=1;
    Ksw1(1)=0.11E-9*V; Ksw1(n)=1;
    Ksw2(1)=0.11E-9*V; Ksw2(n)=1;
end

y0 = [iPlasmids, iPlasmids, cPlasmids, cPlasmids, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, initialBac];

%Simulations
tspan = [0 300];

if iPlasmids==0
    disp('-------------INPUT=0--------------');
    disp('Ksw1')
    disp(Ksw1);
    disp('Ksw2')
    disp(Ksw2);
    %Figures
    figure(1);
    figure(2);
    figure(3);
    figure(4);
    figure(5);
    figure(6);
    figure(7);
    for i=1:length(Ksw1)
        for j=1:length(Ksw2)
            disp(strcat('Ksw1: ', num2str(Ksw1(i)), '<-->Ksw2: ', num2str(Ksw2(j))));
            deqs=odefunc(Ksw1(i), Ksw2(j), Kdt11, Kdt12);
            [t,y]=ode15s(deqs, tspan, y0);
            figure(1);
            plot(t,y(:,Tox1), 'LineWidth', 1.5); hold on;
            figure(2);
            plot(t,y(:,Tox2), 'LineWidth', 1.5); hold on;
            figure(3);
            plot(t,y(:,Tox3), 'LineWidth', 1.5); hold on
            figure(4);
            plot(t,y(:,Bac), 'LineWidth', 1.5); hold on;
            figure(5);
            plot(t,y(:,AT1), 'LineWidth', 1.5); hold on;
            figure(6);
            plot(t,y(:,AT2), 'LineWidth', 1.5); hold on;
            figure(7);
            plot(t,y(:,T1), 'LineWidth', 1.5); hold on;

        end
    end
    figure(1);
    title('Tox1');
    legend (strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(2))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(2))));
    figure(2);
    title('Tox2');
    legend (strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(2))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(2))));
    figure(3);
    title('Tox3');
    legend (strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(2))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(2))));
    figure(4);
    title('Bac');
    legend (strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(2))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(2))));
    figure(5);
    title('AT1');
    legend (strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(2))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(2))));
    figure(6);
    title('AT2');
    legend (strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(2))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(2))));
    figure(7);
    title('T1');
    legend (strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(1)),'--','Ksw2=',num2str(Ksw2(2))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(1))), strcat('Ksw1=',num2str(Ksw1(2)),'--','Ksw2=',num2str(Ksw2(2))));    
else %iPlasmids~=0
    disp('-------------INPUT=1--------------');
    disp('Kdt11');
    disp(Kdt11);
    disp('Kdt12')
    disp(Kdt12);
    %Figures
    figure(1);
    figure(2);
    figure(3);
    figure(4);
    figure(5);
    figure(6);
    for i=1:length(Kdt11)
        for j=1:length(Kdt12)
            disp(strcat('Kdt11: ', num2str(Kdt11(i)), '<-->Kdt12:', num2str(Kdt12(j))));
            deqs=odefunc(Ksw1, Ksw2, Kdt11(i), Kdt12(j));
            [t,y]=ode15s(deqs, tspan, y0);
            figure(1);
            plot(t,y(:,Tox1), 'LineWidth', 1.5); hold on;
            figure(2);
            plot(t,y(:,Tox2), 'LineWidth', 1.5); hold on;
            figure(3);
            plot(t,y(:,Tox3), 'LineWidth', 1.5); hold on;
            figure(4);
            plot(t,y(:,Bac), 'LineWidth', 1.5); hold on;
            figure(5);
            plot(t,y(:,AT1), 'LineWidth', 1.5); hold on;
            figure(6);
            plot(t,y(:,AT2), 'LineWidth', 1.5); hold on;
            figure(7);
            plot(t,y(:,T1), 'LineWidth', 1.5); hold on;
            figure(8);
            plot(t,y(:,Mb2), 'LineWidth', 1.5); hold on;

        end
    end
    figure(1);
    title('Tox1');
    legend (strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(2))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(2))));
    figure(2);
    title('Tox2');
    legend (strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(2))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(2))));
    figure(3);
    title('Tox3');
    legend (strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(2))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(2))));
    figure(4);
    title('Bac');
    legend (strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(2))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(2))));
    figure(5);
    title('AT1');
    legend (strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(2))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(2))));   
    figure(6);
    title('AT2');
    legend (strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(2))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(2))));
    figure(7);
    title('T1');
    legend (strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(2))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(2))));
    figure(8);
    title('Mb2');
    legend (strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(1)),'--','Kdt12=',num2str(Kdt12(2))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(1))), strcat('Kdt11=',num2str(Kdt11(2)),'--','Kdt12=',num2str(Kdt12(2))));

end