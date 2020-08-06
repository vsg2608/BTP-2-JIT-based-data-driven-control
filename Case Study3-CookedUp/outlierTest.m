load ("./data/batch_norm_data.mat");
Data(100,:,4)=[0.8000, 0.2000, 0.400, 0.8000];
qBatch= 5;                      %Query Batch
qTime= 110;                     %Query Time
size_Profile=30;                %Query Profile size

i_qTime=qTime-size_Profile+1;   %Initial query profile time
qProfile= Data(i_qTime:qTime,:,qBatch); %Query Profile
for b=1:1
    tProfile= Data(:,:,b);
    [rProfile,totalCost,iTime]= TWED(qProfile,tProfile,i_qTime);
    sProfiles(:,:,b)= rProfile;
end


hold off
hold on
subplot(1,2,1);
plot(qProfile(:,3),'r');
hold on
plot(sProfiles(:,3,1),'b');
ylim([0 1]);
title('With Outlier')
ylabel('Profile')
legend('qProfile','simialr')

for b=1:1
    tProfile= Data(:,:,b);
    [rProfile,totalCost,iTime]= TWED2(qProfile,tProfile,i_qTime);
    sProfiles(:,:,b)= rProfile;
end

subplot(1,2,2);
plot(qProfile(:,3),'r');
hold on
plot(sProfiles(:,3,1),'b');
ylim([0 1]);
title('Outlier Removal')
ylabel('Profile')
legend('qProfile','Historicl')
