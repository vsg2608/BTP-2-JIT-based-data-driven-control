clear;
load rawData;
qProfile= [3 2 2 4 1 2 4];
tProfile= [3 2 2 2 2 4 1 2 4 4];
global qSize tSize;
%[totalCost,rProfile]= TWED(qProfile, tProfile);
%TWED_profile= rProfile;
%[minDis,rProfile]= movingWindow(qProfile, tProfile);

conversion=y(:,13);
conversion=conversion.'
qProfile= conversion(190:210);
tProfile= conversion(1:130);

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
%}
    qProfile=horzcat(0,qProfile);
    tProfile=horzcat(0,tProfile);
    [x,qSize]= size(qProfile);
    [x,tSize]= size(tProfile);
    A = zeros(tSize,qSize);
    DPenalty=.5;
    for i=2:tSize
        A(i,1)=Inf;
    end
    for i=2:qSize
        A(1,i)=Inf;
    end

    for r= 2:qSize
        for c= 2:r
            i=r-c+2;
            j=c;
            if((j==2) || (j==qSize))
                dPenalty=0;
            else
                dPenalty= DPenalty;
            end
            temp1= A(i-1,j-1)+abs(qProfile(j)-tProfile(i));
            temp2= A(i-1,j)+ abs(qProfile(j)-tProfile(i-1))+dPenalty;
            A(i,j)= min(temp1,temp2);
        end    
    end
    for r= 1:tSize-qSize-1
        for c= 0:qSize-2
            i= qSize+r-c;
            j= c+2;
            if((j==2) || (j==qSize))
                dPenalty=0;
            else
                dPenalty= DPenalty;
            end
            temp1= A(i-1,j-1)+abs(qProfile(j)-tProfile(i));
            temp2= A(i-1,j)+ abs(qProfile(j)-tProfile(i-1))+dPenalty;
            A(i,j)= min(temp1,temp2);
        end
    end

    for r= 2: qSize
        for c= 0: qSize-r
            i= tSize-c;
            j= r+c;
            if((j==2) || (j==qSize))
                dPenalty=0;
            else
                dPenalty= DPenalty;
            end
            temp1= A(i-1,j-1)+abs(qProfile(j)-tProfile(i));
            temp2= A(i-1,j)+ abs(qProfile(j)-tProfile(i-1))+dPenalty;
            A(i,j)= min(temp1,temp2);
        end
    end

    j=qSize;
    i=tSize;
    
    A
    totalCost=A(i,j);
    rProfile=[];
    while(j~=2 & i~=2)
       
        if(A(i-1,j-1)<=A(i-1,j))
            rProfile= horzcat(tProfile(i),rProfile);
            i=i-1;
            j=j-1;
        else
            i=i-1;
        end
    end
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

