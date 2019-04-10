%Similarity measurement and system identification
clear;
load ("./data/batch_norm_data.mat");

Ts=1;
qBatch= 5;                      %Query Batch
qTime= 101;                     %Query Time
size_Profile=30;                %Query Profile size
T_predicts=[];
Y_predicts=[];
Y_actuals=[];
for itr=1:14
    qTime=35+5*itr;
    i_qTime=qTime-size_Profile+1;   %Initial query profile time
    qProfile= Data(i_qTime:qTime,:,qBatch); %Query Profile

    W=[];
    U= Data(i_qTime:qTime,[1,2],qBatch);    %Inputs
    Y= Data(i_qTime:qTime,3,qBatch);        %Outputs

    data = iddata(Y,U,Ts);              %iddata object
    [sys,x0] = ssest(data,1);           %State Space Model

    prediction_time= 10;                 %Time after qPoint to be predicted
    t = 0:Ts:Ts*(size_Profile-1)+Ts*prediction_time;
    uq= Data(i_qTime:qTime+prediction_time,[1,2],qBatch);
    yq= Data(i_qTime:qTime+prediction_time,3,qBatch);
    [y,x] = lsim(sys,uq',t,x0);
    lastPoint= size_Profile +prediction_time;
    T_predicts(itr)= qTime+prediction_time;
    Y_predicts(itr)= y(lastPoint);
    Y_actuals(itr)= yq(lastPoint);
    plot(t,y);
    hold on;
    plot(t,yq);
end

hold off;
scatter(T_predicts,Y_predicts);
hold on;
scatter(T_predicts,Y_actuals);
err = immse(Y_actuals,Y_predicts)