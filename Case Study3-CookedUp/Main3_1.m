%Similarity measurement and system identification
%With TWED
function[Y_predicts, err]= Main3_1(var,P_Time sProfile)
load ("./data/batch_norm_data.mat");

% rng(seed)
% for j=1:4
%     for i=1:10
%       r=randi([1,350],1,1);         
%      Data(r,:,j)=1.25* Data(r,:,j);
%     end
% end

T_predicts=[];
Y_predicts=[];
Y_actuals=[];
qBatch= 5;                      %Query Batch
qTime= 101;                      %Query Time
size_Profile=15;                %Query Profile size
Ts= dT;                          %Delta Time
prediction_time= P_Time;                 %Time after qPoint to be predicted
forward_Profile=prediction_time;
  
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
            %temp= 0.1*ds*ds + 1*dt*dt + 0.1*db*db;
            W(i)=exp(-temp);
        end
        for i=size_Profile+1:size_Profile+forward_Profile
            dt= abs(iTime-i_qTime);
            db= abs(b-qBatch);
            ds=0;
            
            temp= 1*ds*ds + 0.01*dt*dt + 0.1*db*db;
            %temp= 0.1*ds*ds + 1*dt*dt + 0.1*db*db;
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
        cProfile(i,:)=temp/temp2; % 30*4 row column
    end
    
    input=[1,2];
    output=[var+2];
    U= cProfile(:,input);       %Inputs
    Y= cProfile(:,output);      %Outputs

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
    ylabel('output 2')
end
err = immse(Y_actuals,Y_predicts)
end