% function [targetL,targetC,edge,cortexL]=find_act_M1(subjectName,inv_file,edge_file,TimePoints)
function [targetL]=find_target_inArea(inv_file_name,seeds,radius,plotCort)

%-------Direcry settings
% Base_dir='C:\Documents and Settings\Jun Yao\Subject\';
% inv_file_name=[Base_dir,subjectName,'\results_inv\',inv_file];
% edge_file_name=[Base_dir,subjectName,'\results_inv\',edge_file];

%-------read the required info from the files 
[cortexL,Lcount,LNR]=read_Curry_file3(inv_file_name,'LOCATION',0,0);

%--------Select the points within the area---------------
Tind=[];
[n,tmp]=size(seeds);
[m,tmp]=size(cortexL);
for i=1:n
    dis=sqrt( sum( (cortexL-repmat(seeds(i,:),size(cortexL,1),1)).^2,2 )  );
    index=find(dis<radius);
    remain=[];
    if i>1
        for j=1:length(index)
            if isempty(find(Tind==index(j)))
                remain=[remain,index(j)];
            end
        end
    else
        remain=index';
    end
    Tind=[Tind,remain];
end

targetL=cortexL(Tind,:);
% targetC=contrib(Tind,:);
% edge=edge(2:end,:);
if plotCort
	figure (1)
	plot3(cortexL(:,1),cortexL(:,2),cortexL(:,3),'b.');
	hold on
	plot3(targetL(:,1),targetL(:,2),targetL(:,3),'m+');
    plot3(seeds(:,1),seeds(:,2),seeds(:,3),'y.','MarkerSize',4)
    view(90,90)
end