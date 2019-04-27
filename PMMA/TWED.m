
function [rProfile,totalCost,iTime]=TWED(qProfile,tProfile)
    [qSize,xs]= size(qProfile);
    [tSize,xs]= size(tProfile);
    A = zeros(tSize,qSize);
    DPenalty=inf;
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
            temp1= A(i-1,j-1)+ sqrt(sum((qProfile(j,:) - TProfile(i,:)) .^ 2));
            temp2= A(i-1,j)+ sqrt(sum((qProfile(j,:) - TProfile(i-1,:)) .^ 2))+dPenalty;
            if(temp2<temp1)
                TProfile(i,:)=qProfile(j,:);
                temp1= A(i-1,j-1)+sqrt(sum((qProfile(j,:) - TProfile(i,:)) .^ 2))+dPenalty;
            end
            A(i,j)= temp1;
            cost=A(i,j);
        end
        if(minCost>cost)
            minCost=cost;
            minProfile=TProfile;
        end
    end

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
    rProfile=minProfile(i-qSize+1:i,:);
    iTime=i-qSize+1;
end