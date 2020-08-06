%Data Generation
clear;
no_of_batches=5;
[Data,dT]=generate_data(no_of_batches)

%Data normalization
[Data]= normalize(Data);
save ("./data1/batch_norm_data.mat");

ylabels=["Tj_{sp}", "Tau", "X_{m}", "T_{r}"];
for b=1:no_of_batches
    for i=1:4
        subplot(no_of_batches,4,(b-1)*4+i);
        plot(Data(:,i,b));
        xlabel('No. of Samples') 
        ylabel(ylabels(i)) 
    end
end

%function to generate data for bs batches
function [Data,InitialTime]= generate_data(bs)
    for i= 1:bs
        [Data(:,:,i),InitialTime]= PMMA1(0.0014,37.6,85,90,400,i+3);
    end
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