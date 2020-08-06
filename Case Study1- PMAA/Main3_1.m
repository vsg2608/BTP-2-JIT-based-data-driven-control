
%Similarity measurement and system identification
function [Y_predicts,err]=Main3_1(var,P_Time)
load ("./data/batch_norm_data.mat");
T_predicts=[];
Y_predicts=[];
Y_actuals=[];
qBatch= 5;                      %Query Batch
qTime= 101;                      %Query Time
size_Profile=30;                %Query Profile size
Ts= dT;                          %Delta Time
prediction_time= P_Time;                 %Time after qPoint to be predicted
forward_Profile=prediction_time+7;
  
for itr=1:30
    qTime=35+5*itr;
    i_qTime=qTime-size_Profile+1;   %Initial query profile time
    qProfile= Data(i_qTime:qTime,:,qBatch); %Query Profile
    

    cProfile= qProfile;             %Combined Profile
    wProfiles=[];
    sProfiles=[];
    for b=1:qBatch-1
        tProfile= Data(:,:,b);
        [rProfile,totalCost,iTime]= TWED(qProfile,tProfile,i_qTime);
        rProfile= Data(iTime:iTime+size_Profile+forward_Profile,:,b);
        sProfiles(:,:,b)= rProfile;
        W=[];
        for i=1:size_Profile
            dt= abs(iTime-i_qTime);
            db= abs(b-qBatch);
            ds= sqrt(sum((qProfile(i,:) - rProfile(i,:)) .^ 2));
            temp= 1*ds*ds + 0.01*dt*dt + 0.1*db*db;
            W(i)=exp(-temp);
        end
        for i=size_Profile+1:size_Profile+forward_Profile
            dt= abs(iTime-i_qTime);
            db= abs(b-qBatch);
            ds=0;
            
            temp= 1*ds*ds + 0.01*dt*dt + 0.1*db*db;
            W(i)=exp(-temp);
        end
        wProfiles(:,b)= W;
    end

    %Softmax Implementation
    for i=1:size_Profile
        wProfiles(i,:)= wProfiles(i,:)/sum(wProfiles(i,:));
    end

    cProfile= qProfile;             %Combined Profile
    for i=1:size_Profile+forward_Profile
        temp=0;
        temp2=0;
        for b=1:qBatch-1
            w=wProfiles(i,b);
            temp=temp+sProfiles(i,:,b)*w;
            temp2=temp2+w;
        end
        cProfile(i,:)=temp/temp2;
    end
    
    input=[1,2];
    output=[3,4];
    U= cProfile(:,input);       %Inputs
    Y= cProfile(:,output);      %Outputs

   
    data = iddata(Y,U,Ts,'TimeUnit','minutes');
    [sys,x0] = ssest(data,2);
    %[sys,x0] = ssest(data,2,'Ts',data.Ts);
    

    t = 0:Ts:Ts*(size_Profile-1)+Ts*prediction_time;
    uq= Data(i_qTime:qTime+prediction_time,input,qBatch);
    yq= Data(i_qTime:qTime+prediction_time,output,qBatch);
    [y,x] = lsim(sys,uq',t,x0);
    lastPoint= size_Profile +prediction_time;
    T_predicts(itr)= qTime+prediction_time;
    
    compareOut=var;
    Y_predicts(itr)= y(lastPoint,compareOut);
    Y_actuals(itr)= yq(lastPoint,compareOut);
    plot(t,y);
    hold on;
    plot(t,yq);
end
hold off;
plot(T_predicts,Y_predicts,'-o');
hold on;
plot(T_predicts,Y_actuals,'-o');
legend('Predicted','Actual')
xlabel('Time') 
if(compareOut==1)
    ylabel('Conversion')
else
    ylabel('Tr')
end
err = immse(Y_actuals,Y_predicts)
end