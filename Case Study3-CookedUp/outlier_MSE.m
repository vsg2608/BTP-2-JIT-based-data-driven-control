predictionTime=40;
variable=1;
 seed=10;
%   [prediction1,actual,times]=Main1(variable,predictionTime);
[prediction2]=Main2_0(variable,predictionTime,seed);
  [prediction3]=Main2_1(variable,predictionTime,seed);
 [prediction4]=Main3_1(variable,predictionTime,seed);
  [prediction5]=Main3_0(variable,predictionTime,seed);

% hold off;
% plot(times,actual,'-.ok');
% hold on;
% plot(times,prediction2,'--or','LineWidth',2);
% plot(times,prediction3,':ob','LineWidth',2);
%  legend('actual','TWED','TWED2')
% xlabel('Time') 
% if(variable==1)
%     ylabel('Output1')
% else
%     ylabel('Output2')
% end