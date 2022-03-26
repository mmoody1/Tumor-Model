function Np = TUMOR_ODED2(t,N,flag,s,d1,d2,a,a2,b,b2,c1,c2,c3,c4,sigma,ro,k1,k2,k3)  
   %TUMOR_ODE defines the ODEs for the Tumor Model including the drug
   
   % INPUT parameters:
     %s is constant immune cells source rate (#cells/day)
% sigma is the steepness coefficient (#cells/day)
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
    U = N(3); %chemotherapy drug
    X = N(4); %Normal cell population to study toxic effects on normal tissue

  
    Np(1) = s+(ro*I*T)/(sigma+T)-c1*I*T-d1*I - k1*(1-exp(-U));    %immune cells equation
    Np(2) = a*T*(1-b*T)-c2*I*T-c4*X*T-k2*(1-exp(-U)); %tumor cells equation
    Np(3) = -d2*U; %drug dose equation
    Np(4) = a2*X*(1-b2*X)-c4*X*T-k3*(1-exp(-U)); %normal cells equation

end
