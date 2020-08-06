%Data Generation
clear;
no_of_batches=5;
[Data,dT]=generate_data(no_of_batches)

%Data normalization
[Data]= normalize(Data);
save ("./data/batch_norm_data.mat");

ylabels=["Input1", "Input2", "Output1","Output2"];
for b=1:no_of_batches
    for i=1:4
        subplot(no_of_batches,4,(b-1)*4+i);
        plot(Data(:,i,b));
        xlabel('No. of Samples') 
        ylabel(ylabels(i)) 
    end
end

%function to generatedata for bs batches
function [Data,dT]= generate_data(bs)
    for i= 1:bs
        I=i;
        if(i==1)
            I=i+6;
        end
        [Data(:,:,i),dT]= CookedUp(350,I);
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