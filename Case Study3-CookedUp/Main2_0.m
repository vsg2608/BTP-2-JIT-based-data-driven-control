%Similarity measurement and system identification
function[Y_predicts]= Main2_0(var,P_Time)
load ("./data/batch_norm_data.mat");
% Data(100,:,1)=[0.8000, 0.2000, 0.400, 0.8000];
% Data(130,:,4)=[0.8000, 0.2000, 0.400, 0.8000];
% Data(123,:,3)=[0.8000, 0.2000, 0.400, 0.8000];
% Data(150,:,1)=[0.0018    0.1500    0.2345    0.623];
% Data(70,:,2)=[0.0018    0.1500    0.2345    0.623];
% Data(89,:,4)=[0.0018    0.1500    0.2345    0.623];
% Data(85,:,3)=[0.0018    0.1500    0.2345    0.623];

% rng(seed)
% for j=1:4
%     for i=1:10
%      r=randi([1,350],1,1);
%      Data(r,:,j)=1.25* Data(r,:,j);
%     end
% end

qBatch= 5;                      %Query Batch
qTime= 101;                      %Query Time
size_Profile=15;                %Query Profile sizeT_predicts=[];
Y_predicts=[];
Y_actuals=[];
Ts= dT;
prediction_time= P_Time;                 %Time after qPoint to be predicted

for itr=1:30
    qTime=35+5*itr; %It was 25 before just changed the query points. wait havent changed that in other mains.
    i_qTime=qTime-size_Profile+1;   %Initial query profile time
    qProfile= Data(i_qTime:qTime,:,qBatch); %Query Profile
    

    cProfile= qProfile;             %Combined Profile
    wProfiles=[];
    sProfiles=[];
    for b=1:qBatch-1
        tProfile= Data(:,:,b);
        [rProfile,totalCost,iTime]= TWED(qProfile,tProfile,i_qTime);
        sProfiles(:,:,b)= rProfile;
        W=[];
        for i=1:size_Profile
            dt= abs(iTime-i_qTime);
            db= abs(b-qBatch);
            ds= sqrt(sum((qProfile(i,:) - rProfile(i,:)) .^ 2));
            
            temp= 1*ds*ds + 0.01*dt*dt + 0.1*db*db;
            %temp= 1*ds*ds + 1*dt*dt + 0.1*db*db;
            W(i)=exp(-temp);
        end
        wProfiles(:,b)= W;
    end

    %Softmax Implementation
    for i=1:size_Profile
        wProfiles(i,:)= wProfiles(i,:)/sum(wProfiles(i,:));
    end

    cProfile= qProfile;             %Combined Profile
    for i=1:size_Profile
        temp=0;
        temp2=0;
        for b=1:qBatch-1
            w=wProfiles(i,b);
            temp=temp+sProfiles(i,:,b)*w;
            temp2=temp2+w;
        end
        %Addition of query profile point in combined profile
        qWeight= 1;           %Weight corresponding to query point
        temp=temp+qWeight*qProfile(i,:); 
        temp2= temp2+qWeight;
        cProfile(i,:)=temp/temp2;
    end
    
    input=[1,2];
    output=[var+2];
    U= cProfile(1:size_Profile,input);  %Inputs
    Y= cProfile(1:size_Profile,output);      %Outputs

    data = iddata(Y,U,Ts);
    [sys,x0] = ssest(data,1);

    t = 0:Ts:Ts*(size_Profile-1)+Ts*prediction_time;
    uq= Data(i_qTime:qTime+prediction_time,input,qBatch);
    yq= Data(i_qTime:qTime+prediction_time,output,qBatch);
    [y,x] = lsim(sys,uq',t,x0);
    lastPoint= size_Profile +prediction_time;
    T_predicts(itr)= qTime+prediction_time;
    compareOut=1;
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
    ylabel('Output1')
else
    ylabel('Output2')
end
err = immse(Y_actuals,Y_predicts)
end