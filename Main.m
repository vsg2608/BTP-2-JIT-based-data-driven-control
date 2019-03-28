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

data = iddata(Y,U,5);
sys = ssest(data,2);
yp = predict(sys,data,20);
plot(yp);









