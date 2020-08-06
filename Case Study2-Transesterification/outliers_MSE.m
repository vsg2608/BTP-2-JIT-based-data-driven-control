predictionTime=50;
variable=1;
 [prediction1,err]=Main3_1(predictionTime);
 [prediction2,err]=Main2_0(predictionTime);
 [prediction3,err]=Main2_1(predictionTime);

hold off;
plot(times,actual,'-.ok');
hold on;
plot(times,prediction2,'--or','LineWidth',2);
plot(times,prediction3,':ob','LineWidth',2);
 legend('actual','TWED','TWED2')
xlabel('Time') 
if(variable==1)
    ylabel('Output1')
else
    ylabel('Output2')
end