
function[]=Untitled(seed)
a=[1 2 3; 2 6 7; 8 9 2]
rng(seed)
for j=1:2
    for i=1:3
        
        r=randi([1,3],1,1)
        a(r,j)=[8]
        
    end
end
end