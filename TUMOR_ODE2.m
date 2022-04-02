function Np = TUMOR_ODED2(t,N,flag,s,d1,d2,a,a2,b,b2,c1,c2,c3,c4,sigma,ro,k1,k2,k3)  
 %TUMOR_ODE defines the ODEs for the Tumor Model including the drug
   
 % INPUT parameters:
    % s is constant immune cells source rate (#cells/day)
    % sigma is the steepness coefficient (#cells/day)
    % ro is the tumor recruitment rate of immune cells (1/day) 
    % c1 is the tumor deactivation rate of immune cells (1/cell*day)
    % d1 is the natural death rate of immune cells (1/day)
    % d2
    % a is the intrinsic tumor growth rate (1/day)
    % b is the inverse carrying capacity of tumor population (#cells)
    % b2
    % c2 is the immune kill rate of tumor cells (1/cell*day)
    % c3
    % c4
    % k1
    % k2
    % k3
   
 % OUTPUT:
    % N(1) is the immune cell population = I
    % N(2) is the tumor cell population = T
          
 Np=zeros(size(N));
    I=N(1);     % Immune cell population
    T=N(2);     % Tumor cell population
    U = N(3);   % Chemotherapy drug
    X = N(4);   % Normal cell population to study toxic effects on normal tissue

    Np(1) = s+(ro*N(1)*N(2))/(sigma+N(2))-c1*N(1)*N(2)-d1*N(1) - k1*(1-exp(-N(3)));    % Immune cells equation
    if N(2)>0
        Np(2) = a*N(2)*(1-b*N(2))-c2*N(1)*N(2)-c3*X*N(2)-k2*(1-exp(-N(3)));
    else
        Np(2)=0;
    end
     % Tumor cells equation
      if t>1 && mod(t,2)<0.25
        Np(3) = 5/0.25; 
      else
        Np(3) = -d2*N(3);
      end 

          % Drug dose equation
    %Np(4) = a2*X*(1-b2*X)-c4*X*T-k3*(1-exp(-U));                  % Normal cells equation
    Np(4) = 0;
end
