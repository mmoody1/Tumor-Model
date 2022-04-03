% Tumor model

clear all
clc
clf
close all

%input parameter values
s = 0.1181;         % Constant immune cells source rate (#cells/day)
sigma = 20.19;      % Immune threshold rate
ro = 1.131;         % Tumor recruitment rate of immune cells (1/day) 
c1 = 0.00311;       % Tumor deactivation rate of immune cells (1/cell*day)
d1 = 0.3743;        % Natural death rate of immune cells (1/day)
a = 1.636;          % Intrinsic tumor growth rate (1/day)
b = 0.002;          % Inverse carrying capacity of tumor population (#cells)
c2 = 1;             % Immune kill rate of tumor cells (1/cell*day)


% Plotting Nullclines and indicating fixed points

T1 = linspace(212,1000,100);
T2 = linspace(1,10,100);
T3 = linspace(1,1000,100);


%  Nullclines
Nullcline1a = s.*(sigma + T1)./((c1.*T1 + d1).*(sigma + T1)-ro.*T1); %dI/dt = 0
Nullcline1b = s.*(sigma + T2)./((c1.*T2 + d1).*(sigma + T2)-ro.*T2); %dI/dt = 0
Nullcline2= (a.*(1-b.*T3))./c2 ; %dT/dt = 0

hold on

plot(Nullcline1a, T1,'c');
plot(Nullcline1b, T2,'c');
xlim([0 3]); 
 ylim([0 1000]);
 
 plot(Nullcline2, T3,'k--'); 
 set(gca, 'YScale', 'log'); %log scale for tumor axis
 set(gca, 'XScale', 'linear'); %linear scale for immune axis


%  fixed points indicated as 'o'
plot (0.181778,444.444444444444, 'ro'); 
plot(0.763467,266.667, 'bo');
plot(1.60427,8.2976, 'ro');

% % labeling and legend
xlabel('Immune cell population (X10^6)');
ylabel('Tumor Cell population (X10^6)');
lgd = legend('N1a (dI/dt = 0)', 'N1b (dI/dt = 0)', 'N2 (dT/dt = 0)', 'A(0.181, 444.44) Stable', ' B(0.744, 272.73) Unstable', 'C(1.60427,8.2976) Stable','I vs T');
lgd.FontSize = 5;
title('Nullclines and Trajectory cell growth');    

N10 = linspace(0.00001,0.1,2);      % Initial Immune cell population
N20 = linspace(0.1,15,2);         % Initial tumor cell population - tumor burden (x10^6)

tend = 90;         % Simulation length (days)
%sweeping through initial cell populations for immune and tumor
for i = 1:length(N20)
    for j = 1:length(N10)
        N0=[N10(j) N20(i)];

        [t,N] = ode45('TUMOR_ODE',[0 tend],N0,[],s,d1,a,b,c1,c2,sigma,ro); 
         plot(N(:,1), N(:,2)); %plots immune vs tumor
    end
end
hold off

%Population growth plots
figure
subplot(3,1,1)
plot(t,N(:,1));
    xlabel('time')
    ylabel('Immune cells (X10^6)')
xlim([0 90]);

subplot(3,1,2)

plot(t,N(:,2),'m');
    xlabel('time')
    ylabel('tumor cells (X10^6)')
 xlim([0 90]);

 subplot(3,1,3)
plot(N(:,1),N(:,2),'r');
    xlabel('Immune cells (X10^6)')
    ylabel('Tumor cells (X10^6)')

    %Jacobian matrix and eigen values calculation.

%Fixed points values (I,T)
I_eig = 0.181778;
T_eig = 444.44444; ;

%I_eig = 1.60427;
%T_eig = 8.2976;

%I_eig = 0.763467;
%T_eig = 266.667; 

%I_eig = 1.60427
%T_eig = 7.95588


%Jacobian matrix [dI/dI, dI/dT; dT/dI, dT/dT] 
Jac = [1.131*T_eig/(20.19+T_eig)-0.00311*T_eig-0.3743, 1.131*I_eig/(20.19+T_eig) - (20.19+T_eig)^(-2) - 0.00311*I_eig; -T_eig, 1.636-0.014544*T_eig-I_eig];

lamda = eig(Jac); %lamda(1,2)<0 stable, lamda(1,2)>0 unstable

 


