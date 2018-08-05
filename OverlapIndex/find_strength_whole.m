% function [targetL,targetC,edge,cortexL]=find_act_M1(subjectName,inv_file,edge_file,TimePoints)
function [cortexL,targetC]=find_strength_whole(inv_file_name,TimePoints)

%-------Direcry settings
% Base_dir='C:\Documents and Settings\Jun Yao\Subject\';
% inv_file_name=[Base_dir,subjectName,'\results_inv\',inv_file];
% edge_file_name=[Base_dir,subjectName,'\results_inv\',edge_file];

%-------read the required info from the files 
[cortexL,Lcount,LNR]=read_Curry_file3(inv_file_name,'LOCATION',0,0);
contrib=[];

for i=1:TimePoints
    if i==1
        [contrib_tmp,Ccount,CNR]=read_Curry_file4(inv_file_name,'STRENGTH',1,0);
    elseif i==TimePoints
        [contrib_tmp,Ccount,CNR]=read_Curry_file4(inv_file_name,'STRENGTH',1,1);
    else
        [contrib_tmp,Ccount,CNR]=read_Curry_file4(inv_file_name,'STRENGTH',0,1);
    end
    contrib=[contrib,contrib_tmp];
end


targetC=contrib;
