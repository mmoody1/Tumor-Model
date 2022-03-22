function TUMOR_MAIN_1
%tumor model
clear all
clc
clf

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

subplot(5,1,1)
plot(t,N(:,1));
    xlabel('time')
    ylabel('normal cells')

subplot(5,1,2)
plot(t,N(:,2));
    xlabel('time')
    ylabel('tumor cells')

subplot(5,1,3)
plot(t,N(:,3));
    xlabel('time')
    ylabel('immune cells')
legend('Normal Cells', 'Tumor Cells', 'Immune Cells')


 %%plotting nullcines
figure

T = linspace(0,1,20); 
 N = linspace(0,1,20);
 I = linspace(0,1,20);

 %nullsurface of N
 [T1, I1] = meshgrid(T, I);
Nn= 1/b2 - (c4/(r2*b2).*T);
surf(Nn,T1, I1)
hold on

%nullsurface of T
 [N1, I1] = meshgrid(N, I);
Nt = 1/b2 - c2/(r1*b1).*I - (c3/(r1*b1))*N 
surf(N1, Nt, I1)

 hold on
[N1, T1] = meshgrid(N, T);
Ni =   (s.*(alpha + T))./((c1.*T + d1).*(alpha + T) - ro.*T) %nullsurface of I
 surf(N1, T1, Ni)

surf(N, Nt, I)
xlabel('Normal cells')
ylabel('Tumor cells')
zlabel('Immune cells')
 hold off

% axis('equal')
 


