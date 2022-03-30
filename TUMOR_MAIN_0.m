function TUMOR_MAIN
%  opengl software
%tumor model
clear all
clc
clf
%close all
%input parameter values

s = 0.1181; % Constant immune cells source rate (#cells/day)
sigma = 20.19; % immune threshold rate
ro = 1.131; % Tumor recruitment rate of immune cells (1/day) 
c1 = 0.00311; % Tumor deactivation rate of immune cells (1/cell*day)
d1 = 0.3743; %natural death rate of immune cells (1/day)
a = 1.636; %Intrinsic tumor growth rate (1/day)
b = 0.002; %inverse carrying capacity of tumor population (#cells)
c2 = 1; %Immune kill rate of tumor cells (1/cell*day)


N10 = 0.00001;   % Initial Immune cell population
% N20 = 0.25;        % Initial tumor cell population for a small tumor burden (x10^6)
N20= 10;         % Initial tumor cell population for a large tumor burden
tend = 10;       % Simulation length (time)

N0=[N10 N20];


[t,N] = ode45('TUMOR_ODE',[0 tend],N0,[],s,d1,a,b,c1,c2,sigma,ro);

subplot(3,1,1)
plot(t,N(:,1));
    xlabel('time')
    ylabel('Immune cells (X10^6)')

subplot(3,1,2)
plot(t,N(:,2),'m');
    xlabel('time')
    ylabel('tumor cells (X10^6)')
 
subplot(3,1,3)
plot(N(:,1),N(:,2),'r');
    xlabel('Immune cells (X10^6)')
    ylabel('Tumor cells (X10^6)')

 figure %plots Immune and tumor cell growth
 plot(t,N(:,1));
hold on
plot(t,N(:,2),'m');

xlabel('Time (days)')
ylabel('Number of cells (X10^6)')
legend('Immune cells', 'Tumor cells')
ylim([0 215]);

hold off

 figure
%Plotting Nullclines and indicating fixed points

T1 = linspace(0,1000);

%%Nullcline of Immune cells
Nullcline1 = s.*(sigma + T1)./((c1.*T1 + d1).*(sigma + T1)-ro.*T1); %dI/dt = 0
Nullcline2= (a.*(1-b.*T1))./c2 ; %dT/dt = 0
% Quiver
% [x,y] = meshgrid(0:3, 0:1000);
% u = s+(ro.*x.*y)./(sigma+y)-c1.*x.*y-d1.*x;
% v = a.*x.*(1-b.*y)-c2.*x.*y; 
% quiver(x,y,u,v); 
% hold on
plot(Nullcline1, T1,'k');
xlim([0 3]); 
 ylim([0 1000]); %re-scaled
  hold on
 plot(Nullcline2, T1,'k--'); 
 set(gca, 'YScale', 'log');
 set(gca, 'XScale', 'linear');
  hold on
  plot(N(:,1), N(:,2));
%  

 hold on
 
%  fixed points drawn as an x
plot (0.181778,444.444444444444, 'ro'); 
plot(0.763467,266.667, 'bo');
plot(0.97499,202.02, 'bo');
plot(1.58656,14.2431, 'bo');


% % labeling and legend
xlabel('Immune cell population (X10^6)');
ylabel('Tumor Cell population (X10^6)');
lgd = legend('N1 (dI/dt = 0)', 'N2 (dT/dt = 0)','I vs T', 'A(0.181, 444.44) Stable', ' B(0.744, 272.73) Unstable','C(0.97499,202.02) Unstable', 'D(1.58656,14.2431) Unstable');
lgd.FontSize = 5;
title('Nullclines');
%    
hold off

%end
