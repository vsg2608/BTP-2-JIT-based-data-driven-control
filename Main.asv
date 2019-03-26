%{
%Data Generation
clear;
[Data]= generate_data(5);
save ("./data/batch_raw_data.mat");
%Data normalization
clear;
[Data]= normalize();
save ("./data/batch_norm_data.mat");
%}
clear;
load ("./data/batch_norm_data.mat");

qBatch= 5;
qTime= 21;
qPoint= Data(qTime,:,qBatch);
size_Profile=20
i_qTime=qTime-size_Profile+1;
qProfile= Data(i_qTime:qTime,:,qBatch);
X=[]
U=[]
W=[]
Y=[]
for b=1:qBatch-1
    tProfile= Data(:,:,b);
    [rProfile,totalCost,iTime]= TWED(qProfile,tProfile);
    for i= 1:size_Profile
        X=vertcat(X, rProfile(i,:));
        U=vertcat(U, rProfile(i,11));
        Y=vertcat(Y, rProfile(i,13));
    end
end









%function to generatedata for bs batches
function [Data]= generate_data(bs)
    for i= 1:bs
        [~,~,Data(:,:,i)]= PMMA_DataGeneration(i);
    end
end

% function to normalize raw data
function [norm_data]= normalize()
    clear;
    load ("./data/batch_raw_data.mat");
    [ts,xs,bs]= size(Data);

    for i= 1:xs
        minVal = min(min(Data(:,i,:)));
        maxVal = max(max(Data(:,i,:)));
        norm_data(:,i,:) = (Data(:,i,:) - minVal) / ( maxVal - minVal );
    end 
end

