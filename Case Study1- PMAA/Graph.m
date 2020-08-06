predictionTime=40;
variable=1;
[prediction1,actual,times,err1]=Main1(variable,predictionTime);
[prediction2,err2]=Main2_0(variable,predictionTime);
[prediction3,err3]=Main3_1(variable,predictionTime);



hold off;
plot(times,actual,'-.ok');
hold on;
plot(times,prediction1,'-.or','LineWidth',2);
plot(times,prediction2,'--og','LineWidth',2);
plot(times,prediction3,':ob','LineWidth',2);
legend('actual','Method1','Method2','Method3')
xlabel('Time') 
if(variable==1)
    ylabel('Conversion')
else
    ylabel('Temperature')
end