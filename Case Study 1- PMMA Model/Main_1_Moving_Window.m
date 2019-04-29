%Similarity measurement and system identification
%Moving Window Method in this file
clear;
load ("./data/batch_norm_data.mat");

Ts=InitialTime;                 %Step time
qBatch= 5;                      %Query Batch
qTime= 101;                     %Query Time
size_Profile=30;                %Query Profile size
T_predicts=[];
Y_predicts=[];
Y_actuals=[];
for itr=1:30
    qTime=35+5*itr;
    i_qTime=qTime-size_Profile+1;   %Initial query profile time
    qProfile= Data(i_qTime:qTime,:,qBatch); %Query Profile

    W=[];
    input=[1,2];
    output=[4];
    U= Data(i_qTime:qTime,input,qBatch);    %Inputs
    Y= Data(i_qTime:qTime,output,qBatch);    %Outputs

    data = iddata(Y,U,Ts);              %iddata object
    [sys,x0] = ssest(data,1);           %State Space Model

    prediction_time= 10;                 %Time after qPoint to be predicted
    t = 0:Ts:Ts*(size_Profile-1)+Ts*prediction_time;
    uq= Data(i_qTime:qTime+prediction_time,input,qBatch);
    yq= Data(i_qTime:qTime+prediction_time,output,qBatch);
    [y,x] = lsim(sys,uq',t,x0);
    lastPoint= size_Profile +prediction_time;
    T_predicts(itr)= qTime+prediction_time;
    compareOut=1;
    Y_predicts(itr)= y(lastPoint,compareOut);
    Y_actuals(itr)= yq(lastPoint,compareOut);
    plot(t,y(:,compareOut));
    hold on;
    plot(t,yq(:,compareOut));
end

hold off;
scatter(T_predicts,Y_predicts);
hold on;
line(T_predicts,Y_predicts);
scatter(T_predicts,Y_actuals);
line(T_predicts,Y_actuals);
legend('Predicted','Predicted','Actual','Actual')
xlabel('Time') 
ylabel('Conversion')
err = immse(Y_actuals,Y_predicts)