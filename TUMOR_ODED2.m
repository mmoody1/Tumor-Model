function Np = TUMOR_ODED2(t,N,flag,s,d1,d2,a,b,c1,c2,sigma,ro,k1,k2)  
 %TUMOR_ODE defines the ODEs for the Tumor Model including the drug
   
 % INPUT parameters:
    % s is constant immune cells source rate (#cells/day)
    % sigma is the steepness coefficient (#cells/day)
    % ro is the tumor recruitment rate of immune cells (1/day) 
    % c1 is the tumor deactivation rate of immune cells (1/cell*day)
    % d1 is the natural death rate of immune cells (1/day)
    % d2 is the natural decay rate of the drug
    % a is the intrinsic tumor growth rate (1/day)
    % b is the inverse carrying capacity of tumor population (#cells)
    % c2 is the immune kill rate of tumor cells (1/cell*day)   
    % k1 is Drug toxicity to immune cells
    % k2 is Drug toxicity to tumor cells 
     
 % OUTPUT:
    % N(1) is the immune cell population = I
    % N(2) is the tumor cell population = T
    % N(3) is the drug cycle = U     
 Np=zeros(size(N));
    I=N(1);     % Immune cell population
    T=N(2);     % Tumor cell population
    U = N(3);   % Chemotherapy drug
   
    Np(1) = s+(ro*N(1)*N(2))/(sigma+N(2))-c1*N(1)*N(2)-d1*N(1) - k1*(1-exp(-N(3)));    % Immune cells equation
    
    if N(2)>0    
        Np(2) = a*N(2)*(1-b*N(2))-c2*N(1)*N(2)-k2*(1-exp(-N(3)));% Tumor cells equation
    else
        Np(2)=0; %this prevents the drug from killing the tumor cells to negative population (Because that is unrealistic)
    end            
      
    if t>1 && mod(t,2)<0.25  %determines the  number dosing cycles
        Np(3) = 1/0.25; %drug booster equation
      else
        Np(3) = -d2*N(3); %drug decay equation
      end 

         
end
