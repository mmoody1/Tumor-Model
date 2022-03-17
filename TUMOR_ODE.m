function Np = TUMOR_ODE(t,N,flag,s,d1,r1,r2,b1,b2,c1,c2,c3,c4,alpha,ro)
   %TUMOR_ODE defines the ODEs for the Tumor Model
   
   % INPUT parameters:
     % t is time (days)
     % s is the influx of immune cells when tumor cells present
     % d1 is the death rate of immune cells in the absence of tumors
     % r1 is the growth rate associated with tumor cells
     % r2 is the growth rate associated with normal cells
     % b1 is the reciprocol carrying capacities for tumor cells
     % b2 is the reciprocol carrying capacities for normal cells
     % c1 is the competitive coefficient between tumor on immune cells
     % c2 is the competitive coefficient between immune on tumor cells
     % c3 is the competitive coefficient between normal on tumor cells
     % c4 is the competitive coefficient between tumor on normal cells
     % alpha is the immune threshold rate
     % ro is the immune response rate  
   
    % OUTPUT:
     % N(1) is the normal cell population = X
     % N(2) is the tumor cell population = T
     % N(3) is the immune cell population = I
     
     Np=zeros(size(N));
    X=N(1); %normal cell population
    T=N(2); %tumor cell population
    I=N(3); %immune cell population

    Np(1) = r2*X*(1-b2*X)-c4*T*X;  %normal cells
    Np(2) = r1*T*(1-b1*T)-c2*I*T-c3*T*X; %tumor cells 
    Np(3) = s+((ro*I*T)/(alpha+T))-c1*I*T-d1*I; %immune cells
