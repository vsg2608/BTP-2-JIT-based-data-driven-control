%Similarity measurement and system identification
clear;
load ("./data/batch_norm_data.mat");

Ts=5;
qBatch= 5;                      %Query Batch
qTime= 31;                      %Query Time
size_Profile=30;                %Query Profile size
i_qTime=qTime-size_Profile+1;   %Initial query profile time
qProfile= Data(i_qTime:qTime,:,qBatch); %Query Profile

W=[];
U= qProfile(1:size_Profile,[1,2]);  %Inputs
Y= qProfile(1:size_Profile,3);      %Outputs

data = iddata(Y,U,Ts);              %iddata object
[sys,x0] = ssest(data,3);           %State Space Model

prediction_time= 5;                 %Time after qPoint to be predicted
t = 0:5:145+Ts*prediction_time;
uq= Data(i_qTime:qTime+prediction_time,[1,2],qBatch);
yq= Data(i_qTime:qTime+prediction_time,3,qBatch);
[y,x] = lsim(sys,uq',t,x0);
plot(t,y);
hold on;
plot(t,yq);
