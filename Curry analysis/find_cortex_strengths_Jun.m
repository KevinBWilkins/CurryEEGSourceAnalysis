function [contrib,targetL,targetC] = find_cortex_strengths_Jun(cdr_file_name,cortexL,ROI_ind,TimePoints,plotCort)

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

targetL=cortexL(ROI_ind,:);
targetC=contrib(ROI_ind,:);