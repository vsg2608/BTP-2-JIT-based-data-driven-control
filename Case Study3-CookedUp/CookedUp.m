
function[Data,dT]=CookedUp(itr,seed)
    global IN1 IN2
    M=diag([1 1 1 0]);
    x0=zeros(4,1);
    x0(1)=0.1;
    x0(2)=0.1;
    x0(3)=50;
    x0(4)=0.02;
    options=odeset('Mass',M,'RelTol',1e-2,'AbsTol',1e-2);

    v=0;
    IN1=1;
    IN2=1;
    rng(seed,'twister');
    IN1s = normrnd(IN1,0.2,1,1000);     %Random values of Temprature with 
    IN2s = normrnd(IN2,0.1,1,1000);               %Random values of flow rate of monomer
    %k1=normrnd(1,0.01);
    dT=0.2;
    Data=[];
    for i=1:itr
        Ti=(i-1)*dT;
        Tf=(i)*dT;
        IN1=IN1s(i);
        IN2=IN2s(i);
        if(i>50)
            IN2=IN2+2;
        end
        if(i>100)
            IN2=IN2-2.5;
        end
        if(i>70)
            IN1=IN1+0.7;
        end
        if(i>140)
            IN1=IN1+1.4;
        end
        tspan=(Ti:0.1:Tf);
        [t,y]=ode15s(@Nonlinear,tspan,x0,options);
        [m,~] = size(y)
        x0=y(m,:)
        v=v+x0(3)
        data=[IN1,IN2,x0(3),v];
        Data=vertcat(Data,data);
    end


    plot(Data(:,3));

    function dxdt=Nonlinear(t,x)
        dxdt=zeros(4,1);

        %dxdt(1)=(x(1)*x(2))*k1+IN2*sin(x(3));
        dxdt(1)=(x(1)*x(2))+IN2*sin(x(3));
        dxdt(2)=(x(1)/x(3))+x(2)^2;
        dxdt(3)=IN1*sin(x(1)+x(2));
        dxdt(4)=-x(4)+(x(1)^2+x(2)^2);
    end
end