clear;
load ("./data/batch_raw_data.mat");

qProfile= b4_data(21:30,:);
tProfile= b3_data;
[x,qSize]= size(qProfile);
[x,tSize]= size(tProfile);