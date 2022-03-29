%Jacobian matrix and eigen values.

I = 0.181778;
T = 444.444;

%fixed point values.
%I2 = 0.763467, T2= 266.667
%I3 = 0.97499, T3= 202.02

%Jacobian matrix
J = [1.131*T/(20.19+T)-0.00311*T-0.3743, -1.131*I*T*(20.19+T)^(-2)+I*(20.19+T)^(-1)-0.00311*I; -T, 1.636-(2*0.007272*T)-I];

%S = eig(J) S(1,2)<0 stable fixed point, S(1,2)>0 unstable