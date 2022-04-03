opengl software

% Tumor model with drug

clear all
close all
clc
clf
%close all

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
k2 = 1.636;          % Drug toxicity to tumor cells 


N10 = 0.001;      % Initial Immune cell population (*10^6)
N20= 1;           % Initial tumor cell population (tumor burden) (*10^6)
N30 = 1;          %Initial drug given
tend = 100;          % Simulation length (time)

N0=[N10 N20 N30];

opts = odeset('MaxStep',1e-2);

[t,N] = ode23s('TUMOR_ODED2',[0 tend],N0,[opts],s,d1,d2,a,b,c1,c2,sigma,ro,k1,k2);
                                   
subplot(2,1,1)
plot(t,N(:,1)); %plots immune cell growth over time
    xlabel('time')
    ylabel('Immune cells') 

subplot(2,1,2)
plot(t,N(:,2),'m'); %plots tumor growth over time
    xlabel('time')
    ylabel('tumor cells')

figure %drug over time
plot(t, N(:,3)) %plots drug over time
    xlabel('Time')
    ylabel('Drug delivered')


figure 

%Kill rate of Drug

%subplot(2,1,1)

u= linspace(0,10);
Fu=k2.*(1-exp(-u)); %Fu is the per cell kill rate
hold on
plot(u, Fu);
xlabel('Amount of Drug');
ylabel('Per tumor cell kill rate');

yline(k2,'r--','Saturation level')

%the solution for Q1a

 
%  subplot(2,1,2)
% i = linspace(0,10);
% Fu=k1.*(1-exp(-u));
% plot(i, Fu);
% hold on
% yline(k1,'r--','Saturation level')
% % 
% xlabel('Amount of Drug');
% ylabel('Per immune cell kill rate');
