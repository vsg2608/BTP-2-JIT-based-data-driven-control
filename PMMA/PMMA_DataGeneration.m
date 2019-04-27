%% Non Linear Simulation of Polymerisation of MMA Process under isothermal and batch conditions
% Tf = Total Simulation Time
% n  = Number of Time steps
% M_0 = initial Monomer Conversion
% I_0 = initial initiator Conversion 
% [t,y] returns a vector of Time v/s The states at corresponding time
% x[1:16] in MMADae2(t,x) represents [I ; M ; R ; lamda_0 ; lambda_1 ; lambda_2 ; mu_0;
%...... mu_1 ; mu_2 ; tau_1 ; tau_2 ; volume ; conversion ; kt ; kp ; k_td] 

%%
function [ts,Xms,y0s]= PMMA_DataGeneration(seed)
    global E_d E_p E_td k_d0 k0_tm k0_td0 k0_p0 T R_gas A_1 A_2 A_3 A_4 B_1 B_2 B_3 B_4 MW_m R_li R_lm R_vm f k_d k_tc k0_td k_tm kp_0 kt_0 rho_m rho_p Ts     
    %%%%%Parameters for the MMA Dae system
    
    E_d= 128.45E3;                
    E_p= 18.22E3;                
    E_td = 2.937E3;              
    k_d0= 6.3E16;               
    k0_tm = 279.6E6;
    k0_td0= 588E4;
    k0_p0 = 298.2E2;  

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
    
    %%%%%%%%
    %% Main Code for simulation
    I_0= 0.00258;              %Initial moles of initiator
    M_0= 1E6;                  %Initial moles of monomer
    
    M=diag([1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0]);  %% Mass Matrix form
    y0=zeros(1,16);
    
    rng(seed,'twister');
    Tempratures = normrnd(313.15,2,1,1000);     %Random values of Temprature with 
    R_lms= normrnd(1000,200,1,1000);               %Random values of flow rate of monomer
    ts=[];
    Xms=[];
    Temps=[];
    Rlms=[];
    Tis=[];
    y0s=[];
    for i=1:300
        T=Tempratures(i);
        %T=30+273.15;
        if(i>50)    %Step input
            T=T+10;
        end
        if(i>120)
            T=T-20;
        end
        
        R_lm =R_lms(i);
        Y0=y0;
        
        if(R_lm<0)
            R_lm=0;
        end
        Y0(11)=T; % replacaed 11th parameter to temprature as 10th and 11th were exactly same for all times.
        Y0(10)=R_lm; % replacaed 10th parameter to Rlms
        
        InitialTime=1;
        
        if(i==1)
            Ti=0;
            Tf=Ti+InitialTime;               %Time Interval
        end
        
        if(i~=1)
            Ti=InitialTime+(i-2)*InitialTime;
            Tf=Ti+InitialTime;
        end
        
        n=double((Tf-Ti)*2);             %No of parts
        
        A_1= A(1,1)*(T-273.15) + A(1,2);
        A_2 = A(2,1)*(T-273.15) + A(2,2);
        A_3 = A(3,1)*(T-273.15) + A(3,2);
        A_4 = A(4,1)*(T-273.15) + A(4,2);
        B_1 = B(1,1)*(T-273.15) + B(1,2);
        B_2 = B(2,1)*(T-273.15) + B(2,2);
        B_3 = B(3,1)*(T-273.15) + B(3,2);
        B_4 = B(4,1)*(T-273.15) + B(4,2);
        MW_m = 0.10013;                
        R_li = 0;      % R_li and R_lm and R_vm are kept zero for a batch reactor 
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
        
        if(i==1)
            y0(1)= I_0;
            y0(2)= M_0;
            y0(10)= M_0;
            y0(11)= M_0;
            y0(12)= ((M_0*(MW_m))/(rho_m));
            y0(14)= (kt_0)*exp(A_1);
            y0(15)= (kp_0)*exp(B_1);
            y0(16)= k0_td*exp(A_1);
        end
        
        options=odeset('Mass',M );
        tspan=linspace(Ti,Tf,n);
        [t,y] = ode15s(@MMADae2,tspan,y0,options);
        Xm=y(:,13);
        ts=vertcat(ts,t);
        Xms=vertcat(Xms,Xm);
        Temps=vertcat(Temps,T);
        Temps=vertcat(Temps,T);
        Tis=vertcat(Tis,Ti);
        Tis=vertcat(Tis,Tf);
        Rlms= vertcat(Rlms,R_lm);
        Rlms= vertcat(Rlms,R_lm);
        [m,~] = size(y);
        y0=y(m,:);
        Y0(13)=y0(13);
        Y0(12)= MW_m*(y0(6)+y0(9))/(y0(5)+y0(8));
        Y0(14)= MW_m*(y0(5)+y0(8))/(y0(4)+y0(7));
        Y0(15)= Y0(12)/Y0(14);
        y0s=vertcat(y0s,Y0);
    end
    
    subplot(2,1,1);
    plot(ts,Xms), xlabel('Time(min)'), ylabel('Conversion, Xm'), title('Batch ');
    subplot(2,1,2);
    plot(Tis,Temps),xlabel('Time(min)'), ylabel('Temprature');
    s='graphs-BTP2/run0.png';
    s(16)=int2str(seed);
    %saveas(gcf,s);

     

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

end
%%
 