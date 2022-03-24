function TUMOR_MAIN_0
%tumor model
clear all
clc
clf
%  Np(1) = s+((ro*I*T)/(sigma+T))-c1*I*T-d1*I; %immune cells
%     Np(2) = a*T*(1-b*T)-c2*I*T; %tumor cells
s = 0.1181 %
sigma = 20.19; %
ro = 1.131; %
c1 = 0.00311; %
d1 = 0.3743;
a = 1.636;
b = 0.002;
c2 = 1;



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


%  %%plotting nullcines
% figure
% 
% T = linspace(0,1,20); 
%  N = linspace(0,1,20);
%  I = linspace(0,1,20);
% 
%  %nullsurface of N
%  [T1, I1] = meshgrid(T, I);
% Nn= 1/b2 - (c4/(r2*b2).*T);
% surf(Nn,T1, I1)
% hold on
% 
% %nullsurface of T
%  [N1, I1] = meshgrid(N, I);
% Nt = 1/b2 - c2/(r1*b1).*I - (c3/(r1*b1))*N 
% surf(N1, Nt, I1)
% 
%  hold on
% [N1, T1] = meshgrid(N, T);
% Ni =   (s.*(sigma + T))./((c1.*T + d1).*(sigma + T) - ro.*T) %nullsurface of I
%  surf(N1, T1, Ni)
% 
% surf(N, Nt, I)
% xlabel('Normal cells')
% ylabel('Tumor cells')
% zlabel('Immune cells')
%  hold off
% 
% % axis('equal')
%  
% 

