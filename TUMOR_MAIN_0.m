function TUMOR_MAIN
% opengl software
%tumor model
clear all
clc
clf
close all
%input parameter values

s = 0.1181; % Constant immune cells source rate (#cells/day)
sigma = 20.19; % Steepness coefficient (#cells/day)
ro = 1.131; % Tumor recruitment rate of immune cells (1/day) 
c1 = 0.00311; % Tumor deactivation rate of immune cells (1/cell*day)
d1 = 0.3743; %natural death rate of immune cells (1/day)
a = 1.636; %Intrinsic tumor growth rate (1/day)
b = 0.002; %inverse carrying capacity of tumor population (#cells)
c2 = 1; %Immune kill rate of tumor cells (1/cell*day)


N10 = 0.00001;   % Initial Immune cell population
N20 = 1;        % Initial tumor cell population for a small tumor burden (x10^6)
%N20= 10;         % Initial tumor cell population for a large tumor burden
tend = 150;       % Simulation length (time)

N0=[N10 N20];


[t,N] = ode45('TUMOR_ODE',[0 tend],N0,[],s,d1,a,b,c1,c2,sigma,ro);

subplot(3,1,1)
plot(t,N(:,1));
    xlabel('time')
    ylabel('Immune cells')

subplot(3,1,2)
plot(t,N(:,2),'m');
    xlabel('time')
    ylabel('tumor cells')
 
subplot(3,1,3)
plot(N(:,1),N(:,2),'r');
    xlabel('Immune cells')
    ylabel('Tumor cells')

 figure
%Plotting Nullclines and indicating fixed points

T1 = linspace(0,1000);

%%Nullcline of Immune cells
Nullcline1 = s.*(sigma + T1)./((c1.*T1 + d1).*(sigma + T1)-ro.*T1); %dI/dt = 0
Nullcline2= (a.*(1-b.*T1))./c2 ; %dT/dt = 0
plot(Nullcline1, T1,'m');
xlim([0 3]); %set based on figure 7 in paper and intercepts calculated from Nullcline2
ylim([100 1000]); %chopped off unwanted part of that curve so the graph fits better on the screen
  hold on
 plot(Nullcline2, T1,'r'); 
hold on
%fixed points drawn as an x
plot (0.181778,444.444, 'kx'); 
plot(0.743636,272.727, 'kx');
plot(0.97499,202.02, 'kx');

xlabel('Immune cell population (X10^6)');
ylabel('Tumor Cell population (X10^6)');
legend('N1 (dI/dt = 0)', 'N2 (dT/dt = 0)', '(0.181, 444.44) Fixed point', ' (0.744, 272.73) Fixed point','(0.97499,202.02) Fixed point')
title('Immune cells vs Tumor cells')
hold off


end
