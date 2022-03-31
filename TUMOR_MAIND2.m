function TUMOR_MAIND2
 opengl software
%tumor model with drug
clear all
close all
clc
clf
%close all
%input parameter values

s = 0.1181; % Constant immune cells source rate (#cells/day)
sigma = 20.19; % Steepness coefficient (#cells/day)
ro = 1.131; % Tumor recruitment rate of immune cells (1/day) 
c1 = 0.00311; % Tumor deactivation rate of immune cells (1/cell*day)
d1 = 0.3743; %natural death rate of immune cells (1/day)
d2 = 0.5; %rate of decay of drug (1/day)
a = 1.636; %Intrinsic tumor growth rate (1/day)
a2 = 1; %Intrinsic normal cell growth rate (1/day)
b = 0.002; %inverse carrying capacity of tumor population (#cells)
b2 = 1; %inverse carrying capacity for normal cells
c2 = 1; %Immune kill rate of tumor cells (1/cell*day)
c3 = 1; %competitive coefficient of normal cells on tumor cells.
c4 = 1; %competitive coefficient of tumor cells on normal cells
k1 = 0.05; %drug toxicity to immune cells %k values chosen for ideal conditions
k2 = 0.47; %drug toxicity to tumor cells 
k3 = 0.001; %drug toxicity to normal cells %less immunosuppressive drugs

N10 = 0.00001;   % Initial Immune cell population
N20 = 0.25;        % Initial tumor cell population for a small tumor burden (x10^8)
%N20= 10;         % Initial tumor cell population for a large tumor burden
N30 = 0.5;        %Initial drug at tumor site
N40 = 1; %initial normal cell population
tend = 15;       % Simulation length (time)

N0=[N10 N20 N30 N40];


[t,N] = ode23s('TUMOR_ODED2',[0 tend],N0,[],s,d1,d2,a,a2,b,b2,c1,c2,c3,c4,sigma,ro,k1,k2,k3);

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


figure %Kill rate of Drug

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

figure
plot(N(:,3), N(:,1))
 xlabel('Drug delivered')
 ylabel('Tumor cells')
 ylim([0 0.4]);
 
figure 
plot(N(:,3), N(:,4))
xlabel('Drug delivered')
ylabel('Normal cells')


end
