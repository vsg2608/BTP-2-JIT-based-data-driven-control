
predictionTime=20;
variable=1;
[prediction1,actual,times,err1]=Main1(predictionTime);
[prediction2,err2]=Main2_0(predictionTime);
[prediction3,err3]=Main3_1(predictionTime);
err1
err2
err3


hold off;
plot(times,actual,'-.ok');
hold on;

plot(times,prediction1,'-.or','LineWidth',2);
plot(times,prediction2,'--og','LineWidth',2);
plot(times,prediction3,':ob','LineWidth',2);
legend('actual','Method1','Method2','Method3')
xlabel('Time') 
if(variable==1)
    ylabel('Temperature(T_r)')
else
    ylabel('C')
end