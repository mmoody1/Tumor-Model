function Np = TUMOR_ODE(t,N,flag,s,d1,r1,r2,b1,b2,c1,c2,c3,c4,alpha,ro)
    %TUMOR_ODE defines the ODEs for the Tumor Model

     Np=zeros(size(N));
    X=N(1); %normal cell population
    T=N(2); %tumor cell population
    I=N(3); %immune cell population

    Np(1) = r2*X*(1-b2*X)-c4*T*X;  %normal cells
    Np(2) = r1*T*(1-b1*T)-c2*I*T-c3*T*X; %tumor cells 
    Np(3) = s+(ro*I*T/alpha+T)-c1*I*T-d1*I; %immune cells
