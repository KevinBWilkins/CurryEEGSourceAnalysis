%SUBFUNCTIONS--------------------------------------------------------------

function [contrib,targetL,targetC] = find_cortex_strengths_AC(cdr_file_name,cortexL,edge,TimePoints,plotCort)

contrib=[];
fpos=0;

for i=1:TimePoints
    if i==1
        [contrib_tmp,Ccount,CNR]=read_Curry_file4_AC([cdr_file_name(1:end-4),'_strength.cdr'],'STRENGTH',1,0);
        % disp 'done1'
    elseif i==TimePoints
        [contrib_tmp,Ccount,CNR]=read_Curry_file4_AC([cdr_file_name(1:end-4),'_strength.cdr'],'STRENGTH',1,1);
        % disp 'done2'
    else
        [contrib_tmp,Ccount,CNR]=read_Curry_file4_AC([cdr_file_name(1:end-4),'_strength.cdr'],'STRENGTH',0,1);
        % disp 'done3'
    end
    contrib=[contrib,contrib_tmp];
end

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
targetC=contrib(Tind,:);



