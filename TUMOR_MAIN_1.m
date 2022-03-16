function TUMOR_MAIN_1

clear all

s = 0.33;       % Influx of Immune cells when tumor cells present
d1 = 0.2;       % Death rate of immune cells in the absence of tumors
r1 = 1.5;       % Growth rate associated with tumor cells
r2 = 1.0;       % Growth rate associated with normal cells
b1 = 1.0;       % Reciprocol carrying capacities for tumor cells
b2 = 1.0;       % Reciprocol carrying capacities for normal cells
c1 = 1.0;       % Competitive coefficient between tumor on immune cells
c2 = 0.5;       % Competitive coefficient between immune on tumor cells
c3 = 1.0;       % Competitive coefficient between normal on tumor cells
c4 = 1.0;       % Competitive coefficient between tumor on normal cells
alpha = 0.3;    % Immune threshold rate
ro = 0.01;      % Immune response rate 

N10 = 1;        % Initial normal cell population
N20 = 0.5;        % Initial tumor cell population
N30 = 0;        % Initial immune cell population
tend = 50;      % Simulation length (time)

N0 = [N10 N20 N30];

[t, N] = ode45('TUMOR_ODE',[0 tend], N0,[],s,d1,r1,r2,b1,b2,c1,c2,c3,c4,alpha,ro);

subplot(3,1,1)
plot(t,N(:,1));
    xlabel('time')
    ylabel('normal cells')

subplot(3,1,2)
plot(t,N(:,2));
    xlabel('time')
    ylabel('tumor cells')

subplot(3,1,3)
plot(t,N(:,3));
    xlabel('time')
    ylabel('immune cells')