clear;
E_d= 128.45E3;                
E_p= 18.22E3;                
E_td = 2.937E3;              
k_d0= 6.3E16;               
k0_tm = 279.6E6;
k0_td0= 588E4;
k0_p0 = 298.2E2;  
T=50+273.15;                  
A(1,1)= -2.293E-2;
A(1,2)= 7.973E-1;
A(2,1)= -2.379E-1;
A(2,2)= 8.784;
A(3,1)= 1.385E-1;
A(3,2)= -3.609E1;
A(4,1)= -1.745E-3;
A(4,2)= 2.223;
B(1,1)= -4.020E-3;
B(1,2)= 1.246E-1;
B(2,1)= -3.054E-1;
B(2,2)= 9.734;
B(3,1)= 2.706E-1;
B(3,2)= 8.278;
B(4,1)= 2.493E-1;
B(4,2)= -4.459E1;
R_gas = 8.314;                %% universal gas constant J/mol.K
A_1= A(1,1)*(T-273.15) + A(1,2);
A_2 = A(2,1)*(T-273.15) + A(2,2);
A_3 = A(3,1)*(T-273.15) + A(3,2);
A_4 = A(4,1)*(T-273.15) + A(4,2);
B_1 = B(1,1)*(T-273.15) + B(1,2);
B_2 = B(2,1)*(T-273.15) + B(2,2);
B_3 = B(3,1)*(T-273.15) + B(3,2);
B_4 = B(4,1)*(T-273.15) + B(4,2);
MW_m = 0.10013;               
I_0= 0.0258;              
M_0= 500000;                   
R_li = 0.00000;      % R_li and R_lm and R_vm are kept zero for a batch reactor 
R_vm = 0;
f = 0.58;
k_d = k_d0*exp(-E_d/(R_gas*T));       % 1/s
k_tc = 0.0 ;
k0_td = k0_td0*exp(-E_td/(R_gas*T));  % m.cube per mol per secons
k_tm = k0_tm*exp(-E_p/(R_gas*T));     % m.cube per mol per secons
kp_0 = k0_p0*exp(-E_p/(R_gas*T));     % m.cube per mol per secons
kt_0 = k0_td;                         % m.cube per mol per seconds
rho_m = 966.5 - 1.1*(T-273.15);       % Kg/m^3
rho_p = 1200;                         % Kg/m^3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ti=0;
Tf=120;
n=((Tf-Ti)*2);             %No of parts9
u=100000000;
y0=zeros(16,1);
y0(1)= I_0;
y0(2)= M_0;
y0(3)= 0.01;
y0(4)= 0.01;
y0(10)= M_0;
y0(11)= M_0;
y0(12)= ((M_0*(MW_m))/(rho_m));
y0(14)= (kt_0)*exp(A_1);
y0(15)= (kp_0)*exp(B_1);
y0(16)= k0_td*exp(A_1);
[t,y]= MMA_simulation(T,Tf,Ti,n, u, y0(1), y0(2), y0(3), y0(4), y0(5), y0(6), y0(7), y0(8), y0(9), y0(10), y0(11), y0(12), y0(13), y0(14), y0(15), y0(16),M_0,I_0);

MW=[];
Xm=[];
for i=1:200
    MW(i)=MW_m*(y(i,6)+y(i,9))/(y(i,5)+y(i,8));
    Xm(i)=1- y(i,2)/y(i,11);
end
plot(Xm,MW)