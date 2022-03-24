function Np = TUMOR_ODE_0(t,N,flag,s,d1,c1,a,b,c2,sigma,ro)
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
     % N(1) is the immune cell population = I
     % N(2) is the tumor cell population = T
     
     
     Np=zeros(size(N));
    %X=N(1); %normal cell population
    T=N(1); %tumor cell population
    I=N(2); %immune cell population

    %Np(1) = r2*X*(1-b2*X)-c4*T*X;  %normal cells
    Np(1) = s+((ro*I*T)/(sigma+T))-c1*I*T-d1*I; %immune cells
    Np(2) = a*T*(1-b*T)-c2*I*T; %tumor cells 
    