%Similarity measurement and system identification
%Moving Window Method in this file
function [Y_predicts,Y_actuals,T_predicts,err]=Main1(P_Time)
load ("./data/batch_norm_data.mat");

Ts=dT;                          %Step time
qBatch= 5;                      %Query Batch
qTime= 101;                     %Query Time
size_Profile=15;                %Query Profile size
prediction_time= P_Time;                 %Time after qPoint to be predicted
    
T_predicts=[];
Y_predicts=[];
Y_actuals=[];

for itr=1:30
    qTime=35+5*itr;
    i_qTime=qTime-size_Profile+1;   %Initial query profile time
    qProfile= Data(i_qTime:qTime,:,qBatch); %Query Profile

    W=[];
    input=[1,2];
    output=[3];
    U= Data(i_qTime:qTime,input,qBatch);    %Inputs
    Y= Data(i_qTime:qTime,output,qBatch);    %Outputs

    data = iddata(Y,U,Ts);              %iddata object
    [sys,x0] = ssest(data,1);          %State Space Model

    t = 0:Ts:Ts*(size_Profile-1)+Ts*prediction_time;
    uq= Data(i_qTime:qTime+prediction_time,input,qBatch);
    yq= Data(i_qTime:qTime+prediction_time,output,qBatch);
    [y,x] = lsim(sys,uq',t,x0);
    size(y)
    size(x)
    lastPoint= size_Profile +prediction_time;
    T_predicts(itr)= qTime+prediction_time;
    compareOut=1;
    Y_predicts(itr)= y(lastPoint,compareOut); %predicted value from main1
    Y_actuals(itr)= yq(lastPoint,compareOut); %actual value (in the query batch)
%     plot(t,y(:,compareOut));
%     hold on;
%     plot(t,yq(:,compareOut));
end

hold off;
plot(T_predicts,Y_predicts,'-o');
hold on;
plot(T_predicts,Y_actuals,'-o');
legend('Predicted','Actual')
xlabel('Time') 
if(compareOut==1)
    ylabel('Tr')
else
    ylabel('Tr')
end
err = immse(Y_actuals,Y_predicts)
end