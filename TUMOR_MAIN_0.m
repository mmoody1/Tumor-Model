function TUMOR_MAIN
opengl software
%tumor model
clear all
clc
clf
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
N20 = 1;        % Initial tumor cell population for a large tumor
%N20= 1;         % Initial tumor cell population for a small tumor
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
%Nullclines
I1 = linspace(0,5,11);
T1 = linspace(0,500,11);

%Nullcline of Immune cells
Ni = s.*(sigma + T1)./((c1.*T1).*(sigma + T1) + d1.*(sigma+T1)-ro.*T1);
Ti= (a.*(1-b.*T1))./c2 ;
plot(I1,Ni,'m');
hold on
plot(I1,Ti,'r');
xlabel('Immune cell population');
ylabel('Tumor Cell population');
legend('Immune cells', 'Tumor Cells')
title('Immune cells vs Tumor cells')
plot (2.5,0.818, 'rx'); %still working on nailing that fixed point...also one more i havent been able to highlight yet
plot(4.5,0.1636, 'rx');
hold off


end
