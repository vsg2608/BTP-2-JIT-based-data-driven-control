%Data Generation
clear;
[Data]= generate_data(5);
save ("./data/batch_raw_data.mat");
%Data normalization
clear;
[Data]= normalize();
save ("./data/batch_norm_data.mat");


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