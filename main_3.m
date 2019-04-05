%Similarity measurement and system identification
clear;
load ("./data/batch_norm_data.mat");

qBatch= 5;                      %Query Batch
qTime= 31;                      %Query Time
size_Profile=30;                %Query Profile size
i_qTime=qTime-size_Profile+1;   %Initial query profile time
qProfile= Data(i_qTime:qTime,:,qBatch); %Query Profile

cProfile= qProfile;             %Combined Profile
wProfiles=[];
sProfiles=[];
for b=1:qBatch-1
    tProfile= Data(:,:,b);
    [rProfile,totalCost,iTime]= TWED(qProfile,tProfile);
    sProfiles(:,:,b)= rProfile;
    W=[];
    for i=1:size_Profile
        dt= abs(iTime-qTime);
        db= abs(b-qBatch);
        ds= sqrt(sum((qProfile(i,:) - rProfile(i,:)) .^ 2));
        
        W(i,:)=[exp(-ds*ds),exp(-0.01*dt*dt),exp(-0.1*db*db)];
    end
    wProfiles(:,:,b)= W;
end

%Softmax Implementation
for i=1:3
    wProfiles(:,i,:)= wProfiles(:,i,:)/sum(sum(wProfiles(:,i,:)));
end

cProfile= qProfile;             %Combined Profile
for i=1:size_Profile
    temp=0;
    temp2=0;
    for b=1:qBatch-1
        w=wProfiles(i,1,b);
        temp=temp+sProfiles(i,:,b)*w;
        temp2=temp2+w;
    end
    
    cProfile(i,:)=temp/temp2;
end
U= cProfile(1:size_Profile,[1,2]);  %Inputs
Y= cProfile(1:size_Profile,3);      %Outputs



data = iddata(Y,U,5);
plot(data)
[sys,x0] = ssest(data,3);

t = 0:Ts:Ts*(size_Profile-1)+Ts*prediction_time;
uq= Data(i_qTime:qTime+5,[1,2],qBatch);
yq= Data(i_qTime:qTime+5,3,qBatch);
[y,x] = lsim(sys,uq',t,x0);
plot(t,y);
hold on;
plot(t,yq);

