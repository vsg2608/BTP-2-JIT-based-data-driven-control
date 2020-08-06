function[Data, dT]=PMMA(I_0,M_0,Tj,Tr,itr,seed)
rng(seed,'twister')
%Parameters
R=8.314; %KJ/Kmol.K
rho_p=1200; %Kg/m3
rho_j=997;
U=2.55; %KJ/s.m2.K
A=0.0774; %m2
delta_H= -5.78E4; %KJ/kmol
V_j=1.5E-3; %l
V= 2E-3; %l
MW_I=164.21; %kg/kmol
MW_S=106.17;
MW_M=100.12;
MWm=MW_M;
Cp_p=1.47;%KJ/KgK
Cp_m=1.648;
Cp_s=1.70;
Cs= 5.0 ; %kmol/m3
Cp_j= 4.18;
f=0.53;



%solve ode
M2=diag([1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1]);
y0=zeros(27,1);
y0(1)=normrnd(I_0,0.01*I_0); % 0.0014;
y0(1)
y0(2)=M_0; % 3.76;
y0(26)=Tr+273;  %Reactor Temp
y0(9)=7.0E6*exp(-2.6334E4/(R*y0(26)));
y0(10)=3.0233E13*exp(-1.1700E5/(R*y0(26))); %k0p
y0(11)=(y0(7)*MWm)/rho_p;  %phi_p
y0(12)=exp((2.303*(1-y0(11)))/(0.168-8.21E-6*(y0(26)-387)^2+0.03*(y0(26)-y0(11)))); %D
y0(13)=y0(9)/(1+(y0(9)/(y0(12)*y0(10)))); %kp
y0(14)=1.76E9*exp(-1.1704E4/(R*y0(26))); % kt0
y0(15)=1.4540E20*I_0*exp(-1.4584E5/(R*y0(26))); % k0t
y0(16)=y0(14)/(1+(y0(14)/(y0(12)*y0(15)))); %kt
y0(17)=1.58E15*exp(-1.2874E5/(R*y0(26))); % kd
y0(18)=4.661E9*exp(-7.4479E4/(R*y0(26))); % kfm
y0(19)=1.49E9*exp(-6.6197E4/(R*y0(26))); % kfs
% y0(20)=y0(7)/M_0;  %Conversion
y0(21)=(y0(7)+y0(2))*MW_M+Cs*MW_S+I_0*MW_I;                      %rho_mix
y0(22)=(y0(7)*Cp_p+y0(2)*Cp_m+Cs*Cp_s)/(y0(7)+y0(2)+Cs) ;            %Cp
y0(23) = (I_0-y0(1))/I_0; % X_i 
y0(24) = ((2*f*y0(17)*y0(1))/y0(16)).^(1/2); % r0
y0(25) = y0(13)*y0(2)*y0(24); % Rp    
y0(27)=Tr+273;           %Jacket Temp Tj
% yp0=zeros(27,1)
 
options=odeset('Mass',M2,'RelTol',1e-2,'AbsTol',1e-2); % adjust this if Integration fails

dT=1;      %Minutes

Data=[];

%rng(seed,'twister');
Tj_sp=Tj+273;   %Jacket inlet temprature
tau=480;        %Jacket inverse flow rate
Tempratures = normrnd(Tj_sp,2,1,1000);     %Random values of Temprature with 
taus= normrnd(480,10,1,1000);               %Random values of flow rate of monomer


for i=1:itr
    Tj_sp=Tempratures(i);
    tau=tau+1;
    if(i>70)
        Tj_sp= Tj_sp+15;
    end
    if(i>130)
        Tj_sp= Tj_sp-30;
    end
    Ti=(i-1)*dT;
    Tf=Ti+dT;
    tspan=(Ti*60:10:60*Tf); % in seconds
    [t,y]=ode15s(@PMMADae2,tspan,y0,options);
    [m,~] = size(y);
    y0=y(m,:);
    mw = MW_M*((y(m,8) + y(m,5))./(y(m,7)+y(m,4)));
    mn = MW_M*((y(m,7) + y(m,4))./(y(m,6)+y(m,3)));
    PD=mw./mn;
    data=[Tj_sp, tau, y0(20),y0(26)];
    Data=vertcat(Data,data);
end

 
    function[dydt]=PMMADae2(t,y)
    dydt=zeros(27,1);

    dydt(1)=-y(17)*y(1); %I
    dydt(2)=-(y(13)+y(18))*y(2)*y(3); %M
    dydt(3)=2*f*y(17)*y(1)-y(16)*(y(3)^2); %lambda_0
    dydt(4)=2*f*y(17)*y(1)-(y(13)*y(2)*y(3))+(y(18)*y(2)+y(19)*Cs)*(y(3)-y(4))-y(16)*y(3)*y(4);%lambda_1
    dydt(5)=2*f*y(17)*y(1)-((2*y(4))-y(3))*y(13)*y(2)+(y(18)*y(2)+y(19)*Cs)*(y(3)-y(5))-(y(16)*y(3)*y(5));%lambda_2
    dydt(6)=(y(18)*y(2)+y(19)*Cs)*y(3)+0.5*y(16)*y(3)^2; % mu_0
    dydt(7)=(y(18)*y(2)+y(19)*Cs+y(16)*y(3))*y(4);        %mu_1
    dydt(8)=(y(18)*y(2)+y(19)*Cs)*y(5)+y(16)*y(3)*y(5)+y(16)*y(4)^2;  %mu_2
    dydt(9)= -y(9)+7.0E6*exp(-2.6334E4/(R*y(26))); % kp0
    dydt(10)= -y(10)+3.0233E13*exp(-1.1700E5/(R*y(26))); %k0p
    dydt(11)=-y(11)+(y(7)*MWm)/rho_p;  %phi_p
    dydt(12)=-y(12)+exp((2.303*(1-y(11)))/(0.168-8.21E-6*(y(26)-387)^2+0.03*(y(26)-y(11)))); %D
    dydt(13)=-y(13)+y(9)/(1+(y(9)/(y(12)*y(10)))); %kp
    dydt(14)= -y(14)+1.76E9*exp(-1.1704E4/(R*y(26))); % kt0
    dydt(15)= - y(15)+1.4540E20*I_0*exp(-1.4584E5/(R*y(26))); % k0t
    dydt(16)=-y(16)+y(14)/(1+(y(14)/(y(12)*y(15)))); %kt
    dydt(17)= -y(17)+1.58E15*exp(-1.2874E5/(R*y(26))); % kd
    dydt(18)= -y(18)+4.661E9*exp(-7.4479E4/(R*y(26))); % kfm
    dydt(19)= -y(19)+1.49E9*exp(-6.6197E4/(R*y(26))); % kfs
    dydt(20)=-y(20)+(M_0 - y(2))/M_0;  %Conversion
    
    %energy Balance
    dydt(21)=-y(21)+(y(7)+y(2))*MW_M+Cs*MW_S+y(1)*MW_I;                      %rho_mix
    dydt(22)=-y(22)+(y(7)*Cp_p+y(2)*Cp_m+Cs*Cp_s)/(y(7)+y(2)+Cs) ;            %Cp
    dydt(23)= -y(23) + (I_0-y(1))/I_0;  % X_i
    dydt(24)= - y(24) + abs(((2*f*y(17)*y(1))/y(16)).^(1/2));                        % r0
    dydt(25) = -y(25) + y(13)*y(2)*y(24);                                       % Rp
    dydt(26)= (-delta_H*y(25))/(y(22)*y(21)) -U*A*(y(26)-y(27))/(y(22)*V*y(21));  %Reactor T
    dydt(27)=(Tj_sp-y(27))/tau+(U*A*(y(26)-y(27)))/(Cp_j*V_j*rho_j);           %Jacket Temp Tj
        
    end

        
 end