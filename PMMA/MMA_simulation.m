%% Non Linear Simulation of Polymerisation of MMA Process under isothermal and batch conditions
% Tf = Total Simulation Time
% n  = Number of Time steps
% M_0 = initial Monomer Conversion
% I_0 = initial initiator Conversion 
% [t,y] returns a vector of Time v/s The states at corresponding time
% x[1:16] in MMADae2(t,x) represents [I ; M ; R ; lamda_0 ; lambda_1 ; lambda_2 ; mu_0;
%...... mu_1 ; mu_2 ; tau_1 ; tau_2 ; volume ; conversion ; kt ; kp ; k_td] 

%%
function [t,y]= MMA_simulation(Temp,Tf,Ti,n, u, y01, y02, y03, y04, y05, y06, y07, y08, y09, y010, y011, y012, y013, y014, y015, y016,M_0,I_0)
global E_d E_p E_td k_d0 k0_tm k0_td0 k0_p0 T R_gas A_1 A_2 A_3 A_4 B_1 B_2 B_3 B_4 MW_m R_li R_lm R_vm f k_d k_tc k0_td k_tm kp_0 kt_0 rho_m rho_p Ts     
%%%%%Parameters for the MMA Dae system
E_d= 128.45E3;                
E_p= 18.22E3;                
E_td = 2.937E3;              
k_d0= 6.3E16;               
k0_tm = 279.6E6;
k0_td0= 588E4;
k0_p0 = 298.2E2;  
T=Temp;                
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
% I_0= 0.00258;              
% M_0= 1E6;                   
R_li = 0.00000;      % R_li and R_lm and R_vm are kept zero for a batch reactor 
R_lm = u;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation  of MMA %
%Tf=130;
%n=100;
%I_0= 0.00258;              
%M_0= 1E6; 
[t,y]=NLDaeSoln(Tf,Ti,n, u, y01, y02, y03, y04, y05, y06, y07, y08, y09, y010, y011, y012, y013, y014, y015, y016,I_0,M_0);
Xm=1-y(:,2)/y(:,11);
%plot(t,Xm);
%hold on;
%% function handle for Dae simulation in Ode15s %
function Dae2= MMADae2(t,x)
Dae2=zeros(16,1);

Dae2(1) = -k_d*x(1) + R_li;
Dae2(2) = -(x(15) + k_tm)*((x(4)*x(2))/(x(12))) + R_lm - R_vm - x(15)*((x(3)*x(2))/x(12));
Dae2(3) = 2*f*k_d*x(1) - x(15)*((x(3)*x(2))/(x(12)));
Dae2(4) = x(15)*((x(3)*x(2))/(x(12))) - x(14)*(((x(4))^2)/x(12));
Dae2(5) = x(15)*((x(3)*x(2))/(x(12))) + x(15)*((x(2)*x(4))/x(12)) -x(14)*((x(4)*x(5))/x(12)) + k_tm*((x(2)*(x(4) - x(5)))/x(12)) ;
Dae2(6) = x(15)*((x(3)*x(2))/(x(12))) + x(15)*((x(2)*(x(4) + 2*x(5)))/x(12)) - x(14)*((x(4)*x(6))/x(12)) + k_tm*((x(2)*(x(4) - x(6)))/x(12)) ;
Dae2(7) = k_tm*((x(2)*x(4))/x(12)) +(x(16) + 0.5*k_tc)*(((x(4))^2)/x(12));
Dae2(8) = k_tm*((x(2)*x(5))/x(12)) + x(14)*((x(4)*x(5))/x(12));
Dae2(9) = k_tm*((x(2)*x(6))/x(12)) + x(14)*((x(4)*x(6))/x(12)) + k_tc*(((x(5))^2)/x(12));
Dae2(10) =  R_lm - R_vm ;
Dae2(11) =  R_lm ;
Dae2(12) = MW_m*((1/rho_m)-(1/rho_p))*(-(x(15) + k_tm)*((x(4)*x(2))/(x(12))) + R_lm - R_vm - x(15)*((x(3)*x(2))/x(12)));
Dae2(13)= (-1/M_0)*(-(x(15) + k_tm)*((x(4)*x(2))/(x(12))) + R_lm - R_vm - x(15)*((x(3)*x(2))/x(12)));
Dae2(14) =  -x(14)+(kt_0)*exp(A_1 + A_2*(x(13)) + A_3*((x(13))^2) + A_4*((x(13))^3));
Dae2(15) =  -x(15) +(kp_0)*exp(B_1 + B_2*(x(13)) + B_3*((x(13))^2) + B_4*((x(13))^3));
Dae2(16) =  -x(16)+ k0_td*exp(A_1 + A_2*x(13) + A_3*((x(13))^2) + A_4*((x(13))^3));
end

%%function to solve the original non-linear dae from ode15s with mass
%%matrix form
function [t,y]=NLDaeSoln(Tf, Ti,n, u, y01, y02, y03, y04, y05, y06, y07, y08, y09, y010, y011, y012, y013, y014, y015, y016,I_0,M_0)
M=diag([1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0]);  %% Mass Matrix form
y0=zeros(16,1);
y0(1)= y01;
y0(2)= y02;
y0(3)= y03;
y0(4)= y04;
y0(5)= y05;
y0(6)= y06;
y0(7)= y07;
y0(8)= y08;
y0(9)= y09;
y0(10)= y010;
y0(11)= y011;
y0(12)= y012;
y0(13)= y013;
y0(14)= y014;
y0(15)= y015;
y0(16)= y016;
y0(10)= M_0;
y0(11)= M_0;
y0(12)= ((M_0*(MW_m))/(rho_m));
y0(14)= (kt_0)*exp(A_1);
y0(15)= (kp_0)*exp(B_1);
y0(16)= k0_td*exp(A_1);
options=odeset('Mass',M );
tspan=linspace(Ti,Tf,n);
[t,y] = ode15s(@MMADae2,tspan,y0,options);
end
end
%%