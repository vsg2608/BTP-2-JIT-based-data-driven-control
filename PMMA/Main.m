%Similarity measurement and system identification
clear;
load ("./data/batch_norm_data.mat");

qBatch= 5;
qTime= 31;
qPoint= Data(qTime,:,qBatch);
size_Profile=30;
i_qTime=qTime-size_Profile+1;
qProfile= Data(i_qTime:qTime,:,qBatch);
U=[];
W=[];
Y=[];
for b=1:1
    tProfile= Data(:,:,b);
    [rProfile,totalCost,iTime]= TWED(qProfile,tProfile);
    for i= 1:size_Profile
        U=vertcat(U, rProfile(i,[1,2]));
        Y=vertcat(Y, rProfile(i,3));
    end
end
U=[];
Y=[];
for i= 1:size_Profile
    U=vertcat(U, qProfile(i,[1,2]));
    Y=vertcat(Y, qProfile(i,3));
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
