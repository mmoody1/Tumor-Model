% Tumor model with drug
opengl software

clear all
close all
clc
clf

% Input parameter values
s = 0.1181;         % Constant immune cells source rate (#cells/day)
sigma = 20.19;      % Steepness coefficient (#cells/day)
ro = 1.131;         % Tumor recruitment rate of immune cells (1/day) 
c1 = 0.00311;       % Tumor deactivation rate of immune cells (1/cell*day)
d1 = 0.3743;        % Natural death rate of immune cells (1/day)
d2 = 2.0;           % Natural rate of decay of drug (1/day)
a = 1.636;          % Intrinsic tumor growth rate (1/day)
b = 0.002;          % Inverse carrying capacity of tumor population (#cells)
c2 = 1;             % Immune kill rate of tumor cells (1/cell*day)
k1 = 0.05;          % Drug toxicity to immune cells 
k2 = 1.636;         % Drug toxicity to tumor cells 
f = 3;              % dosing frequency
beta = 1;           % booster concentration

% Initial conditions
N10 = 0.001;        % Initial Immune cell population (*10^6)
N20= 1;             % Initial tumor cell population (tumor burden) (*10^6)
N30 = 1;            %Initial drug given
tend = 100;         % Simulation length (time)

N0=[N10 N20 N30];
opts = odeset('MaxStep',1e-2);
[t,N] = ode23s('TUMOR_ODED2',[0 tend],N0,[opts],s,d1,d2,a,b,c1,c2,sigma,ro,k1,k2,f,beta);
                                   
subplot(3,1,1)
plot(t,N(:,1)); %plots immune cell growth over time
    xlabel('time')
    ylabel('Immune cells') 

subplot(3,1,2)
plot(t,N(:,2),'m'); %plots tumor growth over time
    xlabel('time')
    ylabel('Tumor cells')

subplot(3,1,3)
plot(t, N(:,3), 'c') %plots drug over time
    xlabel('Time')
    ylabel('Drug delivered')

% Plot for Amount of Drug vs the tumor cell kill rate k2
figure 
u= linspace(0,10);
Fu=k2.*(1-exp(-u)); %Fu is the per cell kill rate
hold on
plot(u, Fu);
xlabel('Amount of Drug');
ylabel('Per tumor cell kill rate');
yline(k2,'r--','Saturation level')

% Nullcline of Immune and Tumor cells
T1 = linspace(212,1000,100);
T2 = linspace(1,10,100);
T3 = linspace(1,1000,100);

Nullcline1a = s.*(sigma + T1)./((c1.*T1 + d1).*(sigma + T1)-ro.*T1); %dI/dt = 0
Nullcline1b = s.*(sigma + T2)./((c1.*T2 + d1).*(sigma + T2)-ro.*T2); %dI/dt = 0
Nullcline2= (a.*(1-b.*T3))./c2 ; %dT/dt = 0

figure
hold on
plot(Nullcline1a, T1,'c');
plot(Nullcline1b, T2,'c');
xlim([0 3]); 
 ylim([0 1000]); %re-scaled
 plot(Nullcline2, T3,'k--'); 
 set(gca, 'YScale', 'log');
 set(gca, 'XScale', 'linear');
 % % labeling and legend
xlabel('Immune cell population (X10^6)');
ylabel('Tumor Cell population (X10^6)');
lgd = legend('N1 (dI/dt = 0)', 'N2 (dT/dt = 0)','I vs T', 'A(0.181, 444.44) Stable', ' B(0.744, 272.73) Unstable','C(0.97499,202.02) Unstable', 'D(1.58656,14.2431) Unstable');
lgd.FontSize = 5;
title('Nullclines');   

% Trajectory plot lines 
plot(N(:,1),N(:,2), 'm');
lgd = legend('N1 (dI/dt = 0)', 'N1 (dI/dt = 0)','N2 (dT/dt = 0)','I vs T');
lgd.FontSize = 5;
title('Nullclines and Trajectories of treatment');  
hold off
