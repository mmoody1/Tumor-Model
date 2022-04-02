%  opengl software
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

% N
%  ullcline of Immune cells
Nullcline1a = s.*(sigma + T1)./((c1.*T1 + d1).*(sigma + T1)-ro.*T1); %dI/dt = 0
Nullcline1b = s.*(sigma + T2)./((c1.*T2 + d1).*(sigma + T2)-ro.*T2); %dI/dt = 0
Nullcline2= (a.*(1-b.*T3))./c2 ; %dT/dt = 0

hold on

plot(Nullcline1a, T1,'c');
plot(Nullcline1b, T2,'c');
xlim([0 3]); 
 ylim([0 1000]); %re-scaled
  hold on
 plot(Nullcline2, T3,'k--'); 
 set(gca, 'YScale', 'log');
 set(gca, 'XScale', 'linear');
  hold on

% % labeling and legend
xlabel('Immune cell population (X10^6)');
ylabel('Tumor Cell population (X10^6)');
lgd = legend('N1 (dI/dt = 0)', 'N2 (dT/dt = 0)','I vs T', 'A(0.181, 444.44) Stable', ' B(0.744, 272.73) Unstable','C(0.97499,202.02) Unstable', 'D(1.58656,14.2431) Unstable');
lgd.FontSize = 5;
title('Nullclines');    

N10 = linspace(0.1,0.5,3);      % Initial Immune cell population
N20 = linspace(0.1,20,3);         % Initial tumor cell population for a small tumor burden (x10^6)
%N20 = 10;          % Initial tumor cell population for a large tumor burden
tend = 150;         % Simulation length (time)

for i = 1:length(N20)
    for j = 1:length(N10)
        N0=[N10(j) N20(i)];

        [t,N] = ode45('TUMOR_ODE',[0 tend],N0,[],s,d1,a,b,c1,c2,sigma,ro); 
        plot(N(:,1), N(:,2));
    end
end
hold off

