%Similarity measurement and system identification
clear;
load ("./data/batch_norm_data.mat");

qBatch= 5;                      %Query Batch
qTime= 101;                      %Query Time
size_Profile=30;                %Query Profile sizeT_predicts=[];
Y_predicts=[];
Y_actuals=[];
Ts= 1;
for itr=1:30
    qTime=35+5*itr;
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
            dt= abs(iTime-i_qTime);
            db= abs(b-qBatch);
            ds= sqrt(sum((qProfile(i,:) - rProfile(i,:)) .^ 2));
            
            temp= 1*ds*ds + 0.01*dt*dt + 0.1*db*db;
            W(i)=exp(-temp);
        end
        wProfiles(:,b)= W;
    end

    %Softmax Implementation
    for i=1:size_Profile
        wProfiles(i,:)= wProfiles(i,:)/sum(wProfiles(i,:));
    end

    cProfile= qProfile;             %Combined Profile
    for i=1:size_Profile
        temp=0;
        temp2=0;
        for b=1:qBatch-1
            w=wProfiles(i,b);
            temp=temp+sProfiles(i,:,b)*w;
            temp2=temp2+w;
        end
        %Addition of query profile point in combined profile
        qWeight= 1;           %Weight corresponding to query point
        temp=temp+qWeight*qProfile(i,:); 
        temp2= temp2+qWeight;
        cProfile(i,:)=temp/temp2;
    end
    U= cProfile(1:size_Profile,[1,2]);  %Inputs
    Y= cProfile(1:size_Profile,3);      %Outputs

    data = iddata(Y,U,Ts);
    [sys,x0] = ssest(data,1);

    prediction_time= 10;                 %Time after qPoint to be predicted
    t = 0:Ts:Ts*(size_Profile-1)+Ts*prediction_time;
    uq= Data(i_qTime:qTime+prediction_time,[1,2],qBatch);
    yq= Data(i_qTime:qTime+prediction_time,3,qBatch);
    [y,x] = lsim(sys,uq',t,x0);
    lastPoint= size_Profile +prediction_time;
    T_predicts(itr)= qTime+prediction_time;
    Y_predicts(itr)= y(lastPoint);
    Y_actuals(itr)= yq(lastPoint);
    plot(t,y);
    hold on;
    plot(t,yq);
end

hold off;
scatter(T_predicts,Y_predicts);
hold on;
scatter(T_predicts,Y_actuals);
err = immse(Y_actuals,Y_predicts)