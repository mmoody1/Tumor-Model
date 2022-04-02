
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
d2 = 2.0;           % Rate of decay of drug (1/day)
a = 1.636;          % Intrinsic tumor growth rate (1/day)
a2 = 1;             % Intrinsic normal cell growth rate (1/day)
b = 0.002;          % Inverse carrying capacity of tumor population (#cells)
b2 = 1;             % Inverse carrying capacity for normal cells
c2 = 1;             % Immune kill rate of tumor cells (1/cell*day)
c3 = 1;             % Competitive coefficient of normal cells on tumor cells.
c4 = 1;             % Competitive coefficient of tumor cells on normal cells
k1 = 0.05;          % Drug toxicity to immune cells %k values chosen for ideal conditions
k2 = 1.636;          % Drug toxicity to tumor cells 
k3 = 0.001;         % Drug toxicity to normal cells %less immunosuppressive drugs

N10 = 0.001;      % Initial Immune cell population
%N20 = 0.25;         % Initial tumor cell population for a small tumor burden (x10^8)
N20= 1;           % Initial tumor cell population for a large tumor burden
N30 = 5;          %Initial drug at tumor site
N40 = 0;            %initial normal cell population
tend = 100;          % Simulation length (time)

N0=[N10 N20 N30 N40];

opts = odeset('MaxStep',1e-2);

[t,N] = ode23s('TUMOR_ODED2',[0 tend],N0,[opts],s,d1,d2,a,a2,b,b2,c1,c2,c3,c4,sigma,ro,k1,k2,k3);
                                   
subplot(3,1,1)
plot(t,N(:,1)); %plots immune cell growth over time
    xlabel('time')
    ylabel('Immune cells') 

subplot(3,1,2)
plot(t,N(:,2),'m'); %plots tumor growth over time
    xlabel('time')
    ylabel('tumor cells')
    
subplot(3,1,3)
plot(t,N(:,4),'k'); %plots normal cell growth over time
    xlabel('time')
    ylabel('normal cells')



figure %drug over time
plot(t, N(:,3)) %plots drug over time
    xlabel('Time')
    ylabel('Drug delivered')


figure 

%Kill rate of Drug

%subplot(3,1,1)

u= linspace(0,10);
Fu=k2.*(1-exp(-u)); %Fu is the per cell kill rate
plot(u, Fu);
xlabel('Amount of Drug');
ylabel('Per tumor cell kill rate');

hold on
yline(k2,'r--','Saturation level(k = 0.47)')

% subplot(3,1,2)
% 
% Fu=k3.*(1-exp(-u)); %Fu is the per cell kill rate
% plot(u, Fu)
% hold on
% xlabel('Amount of Drug');
% ylabel('Per normal cell kill rate)');
% yline(k3,'r--','Saturation level(k = 0.1)')
% 
% subplot(3,1,3)
% Fu=k1.*(1-exp(-u));
% plot(u, Fu);
% hold on
% yline(k1,'r--','Saturation level(k = 0.2)')
% 
% xlabel('Amount of Drug');
% ylabel('Per immune cell kill rate)');

% figure
% plot(N(:,3), N(:,1))
%  xlabel('Drug delivered')
%  ylabel('Tumor cells')
%  ylim([0 0.4]);
%  
% figure 
% plot(N(:,3), N(:,4))
% xlabel('Drug delivered')
% ylabel('Normal cells')
