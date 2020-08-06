function[Data1]=RandomOutliers(seed)
load ("./data/batch_norm_data.mat");

rng(seed)
for j=1:4
    for i=1:6
     r=randi([1,400],1,1);
     Data(r,:,j)=[0.8120 1 0.0002];
    end
end
end