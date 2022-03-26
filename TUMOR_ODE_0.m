function Np = TUMOR_ODE(t,N,flag,s,d1,a,b,c1,c2,sigma,ro)  
   %TUMOR_ODE defines the ODEs for the Tumor Model
   
   % INPUT parameters:
     %s is constant immune cells source rate (#cells/day)
% alpha is the steepness coefficient (#cells/day)
% ro is the tumor recruitment rate of immune cells (1/day) 
% c1 is the tumor deactivation rate of immune cells (1/cell*day)
% d1 is the natural death rate of immune cells (1/day)
% a is the intrinsic tumor growth rate (1/day)
% b is the inverse carrying capacity of tumor population (#cells)
% c2 is the immune kill rate of tumor cells (1/cell*day)
   
% OUTPUT:
% N(1) is the immune cell population = I
% N(2) is the tumor cell population = T
     
     
Np=zeros(size(N));
    T=N(1); %tumor cell population
    I=N(2); %immune cell population

  
    Np(1) = s+(ro*I*T)/(sigma+T)-c1*I*T-d1*I;    %immune cells equation
    Np(2) = a*T*(1-b*T)-c2*I*T;                    %tumor cells equation

end
