predictionTime=50;
variable=2;
itr=0;
sProfiles=[]
errs1=[]
errs2=[]
errs3=[]
for sProfile=30:35
    [prediction1,actual,times, err1]=Main1(variable,predictionTime, sProfile);
    [prediction2, err2]=Main2_0(variable,predictionTime, sProfile);
    [prediction3, err3]=Main3_1(variable,predictionTime, sProfile);
    sProfiles[itr]= sProfile;
    errs1[itr]= err1;
    errs2[itr]= err2;
    errs3[itr]= err3;
    itr=itr+1;
end

plot(sProfiles,errs1,'-.or','LineWidth',2);
plot(sProfiles,errs2,'--og','LineWidth',2);
plot(sProfiles,errs3,':ob','LineWidth',2);
legend('Main1','Main2','Main3');
xlabel('Profile Size');