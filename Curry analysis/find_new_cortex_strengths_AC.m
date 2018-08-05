function [targetL,targetC] = find_new_cortex_strengths_AC(cortexL,contrib,new_edge,TimePoints,plotCort)

fpos=0;

%--------Select the points inside the edge---------------
Tind=[];
[n,tmp]=size(new_edge);
[m,tmp]=size(cortexL);

for i=1:n
    dis= sqrt(sum((cortexL-repmat(new_edge(i,:),size(cortexL,1),1)).^2,2));
    
    %ind = find(dis<5);
    [yind,ind] = min(dis);
    if dis(ind)<5
        if i==1
            Tind = ind;
        end
        Tind = [Tind;ind];
    end
end
Tind = unique(Tind);

targetL=cortexL(Tind,:);
targetC=contrib(Tind,:);

end