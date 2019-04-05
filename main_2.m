%Similarity measurement and system identification
clear;
load ("./data/batch_norm_data.mat");

qBatch= 5;                      %Query Batch
qTime= 31;                      %Query Time
size_Profile=30;                %Query Profile size
i_qTime=qTime-size_Profile+1;   %Initial query profile time
qProfile= Data(i_qTime:qTime,:,qBatch); %Query Profile

cProfile= qProfile;             %Combined Profile
wProfiles=[]
sProfiles=[]
for b=1:qBatch-1
    tProfile= Data(:,:,b);
    [rProfile,totalCost,iTime]= TWED(qProfile,tProfile);
    sProfiles(:,:,b)= rProfile;
    W=[]
    for i=1:size_Profile
        
    end
end

data = iddata(Y,U,5);
[sys,x0] = ssest(data,3);

t = 0:5:145+5*5;
uq= Data(i_qTime:qTime+5,[1,2],qBatch);
yq= Data(i_qTime:qTime+5,3,qBatch);
[y,x] = lsim(sys,uq',t,x0);
plot(t,y);
hold on;
plot(t,yq);
