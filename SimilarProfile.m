clear;
load rawData;
qProfile= [3 2 2 4 1 2 4];
tProfile= [0 0 3 2 20 4 1 2 4 4];
global qSize tSize;
%[totalCost,rProfile]= TWED(qProfile, tProfile);
%TWED_profile= rProfile;
%[minDis,rProfile]= movingWindow(qProfile, tProfile);

conversion=y(:,13);
conversion=conversion.'
qProfile= conversion(190:210);
tProfile= conversion(170:230);
tProfile(30)=1.0;
[totalCost,rProfile,A]= TWED(qProfile, tProfile);
totalCost
TWED_profile= rProfile;
[minDis,rProfile]= movingWindow(qProfile, tProfile);
hold on
plot(qProfile);
plot(tProfile);
plot(TWED_profile);
minDis
hold off





function [totalCost,rProfile,A]= TWED(qProfile, tProfile)

    qProfile=horzcat(0,qProfile);
    tProfile=horzcat(0,tProfile);
    [x,qSize]= size(qProfile);
    [x,tSize]= size(tProfile);
    A = zeros(tSize,qSize);
    DPenalty=0.01;
    for i=2:tSize
        A(i,1)=0;
    end
    for i=2:qSize
        A(1,i)=Inf;
    end
    
    itr=1;
    minProfile=[];
    minCost=Inf;
    for r= 2:tSize-qSize+2
        TProfile=tProfile;
        for c =0:qSize-2
            i= r+c;
            j= c+2;
            A(i,j)=itr;
            itr=itr+1;
            
            if((j==2) || (j==qSize))
                dPenalty=0;
            else
                dPenalty= DPenalty;
            end
            temp1= A(i-1,j-1)+abs(qProfile(j)-TProfile(i));
            temp2= A(i-1,j)+ abs(qProfile(j)-TProfile(i-1))+dPenalty;
            if(temp2<temp1)
                TProfile(i)=qProfile(j);
                temp1= A(i-1,j-1)+abs(qProfile(j)-TProfile(i))+dPenalty;
            end
            A(i,j)= temp1;
            cost=A(i,j);
        end
        if(minCost>cost)
            minCost=cost;
            minProfile=TProfile;
        end
    end
    
    minCost;
    minProfile;
    j=qSize;
    i=tSize;
    minCost=Inf;
    for c= qSize:tSize
        if(minCost>A(c,qSize))
            minCost=A(c,qSize);
            i=c;
        end
    end
    totalCost=minCost;
    rProfile=[];
    while(j~=1 & i~=1)
       
        if(A(i-1,j-1)<=A(i-1,j))
            rProfile= horzcat(minProfile(i),rProfile);
            i=i-1;
            j=j-1;
        else
            i=i-1;
        end
    end
    rProfile;
end

%Moving window method
function [minDis, rProfile]= movingWindow(qProfile,tProfile)
    [x,qSize]= size(qProfile);
    [x,tSize]= size(tProfile);
    minDis= Inf;
    index=1;
    for i= 0:(tSize-qSize)
        dist= 0;
        for j = 1:qSize
            dist= dist+ abs(qProfile(j)-tProfile(i+j));
        end
        disp(dist);
        if(minDis>dist)
            minDis=dist;
            index=i+1;
        end
    end
    rProfile= tProfile(index:index+qSize-1);
end

%TWED

