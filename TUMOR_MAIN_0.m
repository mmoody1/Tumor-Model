function TUMOR_MAIN_0
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
N20 = 10;        % Initial tumor cell population for a large tumor
%N20= 1;         % Initial tumor cell population for a small tumor
tend = 50;       % Simulation length (time)

N0 = [N10 N20];

 [t, N] = ode23s('TUMOR_ODE_0',[0 tend], N0,[],s,d1,c1,a,b,c2,sigma,ro);

subplot(2,1,1)
plot(t,N(:,1));
    xlabel('time')
    ylabel('Immune cells')

subplot(2,1,2)
plot(t,N(:,2));
    xlabel('time')
    ylabel('tumor cells')

legend('Immune cells', 'Tumor Cells')

