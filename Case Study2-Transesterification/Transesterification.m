%[t,y] returns a vector of Time v/s The states at corresponding time
%x[1:6] in TRES(t,x) represents [CTG ; CDG ; CMG ; CE ; CA ;CGL ;Tr; Qj; Tc; Fc; r; k1; k2; k3; k4; k5; k6] 
function[Data,dT]=Transesterification(Ts,Tr,itr,seed)
    rng(seed,'twister');
    global a b rho_r V C_w rho_w UA Delta_H V_w Tsp C_r MW_r FdT
    a(1)=3.92E7/60;
    %a(1)=normrnd(3.92E7/60,0.15*3.92E7/60);
    a(2)=5.77E5/60;
    %a(2)=normrnd(5.77E5/60,0.15*5.77E5/60);
    a(3)=5.88E12/60 ;
    %a(3)=normrnd(5.88E12/60,0.15*5.88E12/60); 
    a(4)=0.98E10/60 ;
    %a(4)=normrnd(0.98E10/60,0.15*0.98E10/60);
    a(5)=5.35E3/60 ;
    %a(5)=normrnd(5.35E3/60,0.15*5.35E3/60);
    a(6)=2.15E4/60;
    %a(6)=normrnd(2.15E4/60,0.15*2.15E4/60);
    b(1)=6614.83;
    b(2)=4997.98;
    b(3)=9993.96;
    b(4)=7366.64;
    b(5)=3231.18;
    b(6)=4824.87;

    rho_r= 860; % Kg per metre cube
    V= 1;         % metre cube
    Delta_H= -18500; % KJ/Kmol
    UA=7.5; % KJ per sec per Kelvin
    Tsp= 300; % K
    rho_w= 1000; % Kg per metre cube
    C_w= 4.12 ; %KJ per Kg per Kelvin
    C_r = 1277; %KJ per Kmol per Kelvin
    MW_r= 391.40; % Kg per Kmol
    V_w= 0.5; % metre cube

    % solve the original ode from ode15s
    y0=zeros(17,1);
    y0(1)=0.3226;
    y0(5)= 1.9356;
    y0(7)=Tr; %Initial value of Reactor temp.
    %y0(8)=202.5;
    y0(9)=Tr;
    y0(8)= UA*(y0(7)- y0(9));
    y0(10)=0.0058;          %Kg/s (Flow Rate of Coolant Jacket)
    %y0(11)= 0.0578;         %mol/l.s   
    y0(12)=a(1)*exp(-b(1)/y0(7));
    y0(13)=a(2)*exp(-b(2)/y0(7));
    y0(14)=a(3)*exp(-b(3)/y0(7));
    y0(15)=a(4)*exp(-b(4)/y0(7));
    y0(16)=a(5)*exp(-b(5)/y0(7));
    y0(17)=a(6)*exp(-b(6)/y0(7));
    y0(11)= ((y0(12)*y0(1)*y0(5))-(y0(13)*y0(2)*y0(4))+(y0(14)*y0(2)*y0(5))-(y0(15)*y0(3)*y0(4))+(y0(16)*y0(3)*y0(5))-(y0(17)*y0(6)*y0(4)));
    M=diag([ 1 1 1 1 1 1 1 0 1 1 0 0 0 0 0 0 0]);
    options=odeset('Mass',M);

    dT=10;    
    Data=[];
    %rng(seed,'twister');
    Tsp=Ts;   %Jacket inlet temprature
    FdT=0;      %Rate of change of flow rate
    Tempratures = normrnd(Tsp,2,1,1000);     %Random values of Temprature with 
    FdTs= normrnd(FdT,0.0005,1,1000);               %Random values of flow rate of monomer
    for i=1:itr
        Tsp= Tempratures(i);
        FdT=0;
        if(i==250)
            FdT=FdT+0.005;
        end
        if(i>100)
            Tsp=Tsp+10;
        end
        if(i>150)
            Tsp=Tsp-20;
        end
        Ti=(i-1)*dT;
        Tf=Ti+dT;
        tspan=(Ti:1:Tf); % in seconds
        [t,y]=ode15s(@TRES2,tspan,y0,options);
        [m,~] = size(y);
        y0=y(m,:);
        data=[Tsp,y0(10),y0(7)];
       
        Data=vertcat(Data,data);
    end
    
    %function handle for ODE simulation using ode15
    function dxdt= TRES2(t,x)
    dxdt=zeros(17,1);

    dxdt(1)=(-x(12)*x(1)*x(5))+(x(13)*x(2)*x(4));
    dxdt(2)=(x(12)*x(1)*x(5))-(x(13)*x(2)*x(4))-(x(14)*x(2)*x(5))+(x(15)*x(3)*x(4));
    dxdt(3)=(x(14)*x(2)*x(5))-(x(15)*x(3)*x(4))-(x(16)*x(3)*x(5))+(x(17)*x(6)*x(4));
    dxdt(4)=((x(12)*x(1)*x(5))-(x(13)*x(2)*x(4))+(x(14)*x(2)*x(5))-(x(15)*x(3)*x(4))+(x(16)*x(3)*x(5))-(x(17)*x(6)*x(4)));
    dxdt(5)=-((x(12)*x(1)*x(5))-(x(13)*x(2)*x(4))+(x(14)*x(2)*x(5))-(x(15)*x(3)*x(4))+(x(16)*x(3)*x(5))-(x(17)*x(6)*x(4)));
    dxdt(6)=(x(16)*x(3)*x(5))-(x(17)*x(6)*x(4));
    dxdt(7)= (MW_r/(rho_r*V*C_r))*(-x(11)*V*Delta_H - x(8));
    dxdt(8)= -x(8)+ UA*(x(7)- x(9));
    dxdt(9)= (x(10)/(V_w*rho_w))*(Tsp - x(9)) + (x(8)/(rho_w*C_w*V_w)); 
    dxdt(10)= FdT; %flow rate is taken constant
    dxdt(11) = -x(11) + ((x(12)*x(1)*x(5))-(x(13)*x(2)*x(4))+(x(14)*x(2)*x(5))-(x(15)*x(3)*x(4))+(x(16)*x(3)*x(5))-(x(17)*x(6)*x(4)));
    dxdt(12) = -x(12) + a(1)*exp(-b(1)/x(7));
    dxdt(13) =-x(13) + a(2)*exp(-b(2)/x(7));
    dxdt(14) = -x(14) +a(3)*exp(-b(3)/x(7));
    dxdt(15) = -x(15) + a(4)*exp(-b(4)/x(7));
    dxdt(16)= -x(16) + a(5)*exp(-b(5)/x(7));
    dxdt(17) =-x(17) + a(6)*exp(-b(6)/x(7));
    end
end