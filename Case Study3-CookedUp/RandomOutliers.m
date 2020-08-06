function[Data1]=RandomOutliers(seed)
load ("./data/batch_norm_data.mat");

rng(seed)
for j=1:4
    for i=1:6
     r=randi([1,350],1,1);
     Data(r,:,j)=[0.8,0.1,0.12, 0.001];
    end
end
end
