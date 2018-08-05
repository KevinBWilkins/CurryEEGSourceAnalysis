%SUBFUNCTIONS--------------------------------------------------------------

function [targetL] = find_cortex_locations_AC(cortexL,edge)

%--------Select the points inside the edge---------------
Tind=[];
[n,tmp]=size(edge);
[m,tmp]=size(cortexL);

for i=1:n
    dis= sqrt(sum((cortexL-repmat(edge(i,:),size(cortexL,1),1)).^2,2));
    
%     ind = find(dis<2);      %if want to look around pts
    [yind,ind] = min(dis); %if want to pick exact pt that corresponds to selected pt in Curry // This should be just looking at the minimum distance --jun
    
    if dis(ind)<6
        if i==1
            Tind = ind;
        end
        Tind = [Tind;ind];
    end
end
Tind = unique(Tind);

targetL=cortexL(Tind,:);




