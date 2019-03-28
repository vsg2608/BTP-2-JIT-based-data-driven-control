%Similarity measurement and system identification
clear;
load ("./data/batch_norm_data.mat");

qBatch= 5;
qTime= 21;
qPoint= Data(qTime,:,qBatch);
size_Profile=20;
i_qTime=qTime-size_Profile+1;
qProfile= Data(i_qTime:qTime,:,qBatch);
X=[];
U=[];
W=[];
Y=[];
for b=1:1
    tProfile= Data(:,:,b);
    [rProfile,totalCost,iTime]= TWED(qProfile,tProfile);
    for i= 1:size_Profile
        X=vertcat(X, rProfile(i,:));
        U=vertcat(U, rProfile(i,11));
        Y=vertcat(Y, rProfile(i,13));
    end
end











