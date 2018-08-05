% function [targetL,targetC,edge,cortexL]=find_act_M1(subjectName,inv_file,edge_file,TimePoints)
function [cortexL,targetL,targetC]=find_target_PC2(inv_file_name,edge1_file_name,edge2_file_name,edge3_file_name,edge4_file_name,TimePoints,plotCort)

%-------Direcry settings
% Base_dir='C:\Documents and Settings\Jun Yao\Subject\';
% inv_file_name=[Base_dir,subjectName,'\results_inv\',inv_file];
% edge_file_name=[Base_dir,subjectName,'\results_inv\',edge_file];

%-------read the required info from the files 
[cortexL,Lcount,LNR]=read_Curry_file3(inv_file_name,'LOCATION',0,0);
[edge1,Ecount,ENR]=read_Curry_file3(edge1_file_name,'LOCATION',0,0);
[edge2,Ecount,ENR]=read_Curry_file3(edge2_file_name,'LOCATION',0,0);
[edge3,Ecount,ENR]=read_Curry_file3(edge3_file_name,'LOCATION',0,0);
[edge4,Ecount,ENR]=read_Curry_file3(edge4_file_name,'LOCATION',0,0);
edge=[edge1;edge2;edge3;edge4];
% TimePoints=read_Curry_file3(inv_file,'Time',0,0);
contrib=[];
fpos=0;
% now_read=1;
% timePoint=1;
% timenow=1;
for i=1:TimePoints
    if i==1
        [contrib_tmp,Ccount,CNR]=read_Curry_file4(inv_file_name,'FIELDS',1,0);
    elseif i==TimePoints
        [contrib_tmp,Ccount,CNR]=read_Curry_file4(inv_file_name,'FIELDS',1,1);
    else
        [contrib_tmp,Ccount,CNR]=read_Curry_file4(inv_file_name,'FIELDS',0,1);
    end
    contrib=[contrib,contrib_tmp];
end

%--------Select the points inside the edge---------------
Tind=[];
% center=edge(1,:);
[n,tmp]=size(edge);
[m,tmp]=size(cortexL);
% MVe=edge(2:end,:)-repmat(center,n-1,1);
% for i=1:m
%     Vp=cortexL(i,:)-center;
%     for j=1:n-1
%         costheta(j)=dot(Vp,MVe(j,:))/norm(Vp)/norm(MVe(j,:));
%     end
% %     if ~isempty(find(costheta>0.98))
% %         index=find(costheta==1);
% %     end
% %     [minTheta,index]=max(abs(costheta));
%     [minTheta,index]=max(costheta);
% %     closeMVe=
%     if norm(Vp)<=norm(MVe(index,:))
%         Tind=[Tind,i];
%     end
% end
% cur_L=cortexL;
for i=1:n
    dis=sqrt( sum( (cortexL-repmat(edge(i,:),size(cortexL,1),1)).^2,2 )  );
    index=find(dis<5);
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
%     KpInd=find(dis>=5);

%     for j=1:length(index)
    Tind=[Tind,remain];
%     tmp=cur_L(KpInd,:);
%     clear cur_L
%     cur_L=tmp;
%     clear tmp
end

targetL=cortexL(Tind,:);
targetC=contrib(Tind,:);
% edge=edge(2:end,:);
if plotCort
	figure (1)
	plot3(cortexL(:,1),cortexL(:,2),cortexL(:,3),'b.');
	hold on
	plot3(targetL(:,1),targetL(:,2),targetL(:,3),'m+');
end