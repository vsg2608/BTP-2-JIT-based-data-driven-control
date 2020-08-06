
%Data Generation
clear;
no_of_batches=5;
[Data,dT]=generate_data(no_of_batches);

%Data normalization
[Data]= normalize(Data);
save ("./data/batch_norm_data.mat");

ylabels=["Tj_{sp}", "F_j", "T_{r}"];
for b=1:no_of_batches
    for i=1:3
        subplot(no_of_batches,3,(b-1)*3+i);
        plot(Data(:,i,b));
        xlabel('No. of Samples') 
        ylabel(ylabels(i)) 
    end
end

%function to generatedata for bs batches
function [Data,dT]= generate_data(bs)
    for i= 1:bs
        [Data(:,:,i),dT]= Transesterification(380,400,400,i+1);
    end
end

%function to normalize raw data
function [norm_data]= normalize(Data)
    [ts,xs,bs]= size(Data);
    
    for i= 1:xs
        minVal = min(min(Data(:,i,:)));
        maxVal = max(max(Data(:,i,:)));
        norm_data(:,i,:) = (Data(:,i,:) - minVal) / ( maxVal - minVal );
    end 
end