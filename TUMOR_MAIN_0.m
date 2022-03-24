function TUMOR_MAIN_0
%tumor model
clear all
clc
clf
%  Np(1) = s+((ro*I*T)/(sigma+T))-c1*I*T-d1*I; %immune cells
%     Np(2) = a*T*(1-b*T)-c2*I*T; %tumor cells
s = 0.1181 % Constant immune cells source rate (#cells/day)
sigma = 20.19; % Steepness coefficient (#cells/day)
ro = 1.131; % Tumor recruitment rate of immune cells (1/day) 
c1 = 0.00311; % Tumor deactivation rate of immune cells (1/cell*day)
d1 = 0.3743; %natural death rate of immune cells (1/day)
a = 1.636; %Intrinsic tumor growth rate (1/day)
b = 0.002; %inverse carrying capacity of tumor population (#cells)
c2 = 1; %Immune kill rate of tumor cells (1/cell*day)



N10 = 0.5;        % Initial Immune cell population
N20 = 10;        % Initial tumor cell population
tend = 50;      % Simulation length (time)

N0 = [N10 N20];

[t, N] = ode45('TUMOR_ODE_0',[0 tend], N0,[],s,ro,sigma,c1,d1,a,b,c2);

subplot(2,1,1)
plot(t,N(:,1));
    xlabel('time')
    ylabel('Immune cells')

subplot(2,1,2)
plot(t,N(:,2));
    xlabel('time')
    ylabel('tumor cells')

%legend('Immune cells', 'Tumor Cells')


