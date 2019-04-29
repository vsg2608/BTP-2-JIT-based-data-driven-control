%Data Generation
clear;
no_of_batches=5;
[Data,InitialTime]= generate_data(no_of_batches);
save ("./data/batch_raw_data.mat");
 
%Data ploting

ylabels=["Flow Rate", "Temprature", "Avg Mol wt", "Conversion"];
for b=1:no_of_batches
    for i=1:4
        subplot(no_of_batches,4,(b-1)*4+i);
        plot(Data(:,i,b));
        xlabel('Time') 
        ylabel(ylabels(i)) 
    end
end

%Data normalization
[Data]= normalize(Data);
save ("./data/batch_norm_data.mat");

%function to generatedata for bs batches
function [Data,InitialTime]= generate_data(bs)
    for i= 1:bs
        [~,~,Data(:,:,i),InitialTime]= PMMA_DataGeneration(i+1);
    end
    Data= Data(:,[10,11,12,13,14,15],:); % Only 10th(Rlm), 11th(Temp) and 13th(conversion) required.
end

% function to normalize raw data
function [norm_data]= normalize(Data)
    [ts,xs,bs]= size(Data);

    for i= 1:xs
        minVal = min(min(Data(:,i,:)));
        maxVal = max(max(Data(:,i,:)));
        norm_data(:,i,:) = (Data(:,i,:) - minVal) / ( maxVal - minVal );
    end 
end