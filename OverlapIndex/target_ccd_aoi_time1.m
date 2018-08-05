% This function read the CCD information in the target local area.
function ActA=target_ccd_aoi_time1(subject)
% function target_ccd (subject,cdr_file,edge_file,timePoints,plotCort)
% function plot_target_ccd (subjectName,cdr_file,edge_file,timePoints,plotCort)
method={'1','2.5','3','4','5','6','7'};
pre_method=method(1);
addpath('C:\Documents and Settings\Jun Yao\Matlab\Eeg\inv')

% subject={'brad';'Jules';'TK';'Ty';'LK'};
% subject={'BH';'JD';'Ty';'LK'};
% subject={'JD'};

Base_dir='C:\Documents and Settings\Jun Yao\Subject\';
subNo=length(subject);
% taskList={'ext';'flex';'sabd';'sadd'};
% taskList={'ext';'flex';'abd';'add'};
taskList={'abd';'flex'};
% taskList={'ind';'sld'};
% method='LRT15';
% phaseList=[-700, -250, 600;  -700, -250, 600; -700, -250, 600; -700, -150, 700; -700, -200, 600];
%phaseList={'pre'; 'mov'; 'hd'};  
phaseList={'whole'};

% phaseList=[ -700, -200, 600];
plotCort=0;

show_frame=0;
yjcolor=['b','r','m'];
currents={};
for k=1:subNo
    subjectName=subject{k};
	edge1_file_name=[Base_dir,subjectName,'\Data\MRI\','SM1_L.sp'];
    edge2_file_name=[Base_dir,subjectName,'\Data\MRI\','PREMOTOR_L.sp'];
	edge3_file_name=[Base_dir,subjectName,'\Data\MRI\','SM1_R.sp'];
    edge4_file_name=[Base_dir,subjectName,'\Data\MRI\','PREMOTOR_R.sp'];
	
    [edge1,Ecount,ENR]=read_Curry_file3(edge1_file_name,'LOCATION',0,0);
	[edge2,Ecount,ENR]=read_Curry_file3(edge2_file_name,'LOCATION',0,0);
	[edge3,Ecount,ENR]=read_Curry_file3(edge3_file_name,'LOCATION',0,0);
	[edge4,Ecount,ENR]=read_Curry_file3(edge4_file_name,'LOCATION',0,0);
	edge=[edge1;edge2;edge3;edge4];
    
    censul_file_name=[Base_dir,subjectName,'\Data\MRI\','CS_L.sp'];
    [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name,'LOCATION',0,0);
    censul_file_name=[Base_dir,subjectName,'\Data\MRI\','CS_R.sp'];
    [censul2,Ecount,ENR]=read_Curry_file3(censul_file_name,'LOCATION',0,0);
    censul=[censul1;censul2];
    record=[];
    for i=1:length(taskList)
        cur_task=taskList{i};
        for phTime=1:length(phaseList)
            timePoints=13;
            cdr_file=[cur_task,'_',method{phTime},'.cdr'];
%             cdr_file=[cur_task,'.cdr'];

            disp('********************');
            disp(cdr_file)
			cdr_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_file];
			
% 			[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
			[cortexL,targetL,targetC]=find_target_PC3(cdr_file_name,edge,timePoints,plotCort);

            minX=min(targetL(:,1));
            maxX=max(targetL(:,1));
            minY=min(targetL(:,2));
            maxY=max(targetL(:,2));
            minZ=min(targetL(:,3));
            maxZ=max(targetL(:,3));
            x=[minX-5:0.8:maxX+5];
            y=[minY-5:0.8:maxY+5];
            [XI,YI] = meshgrid(x,y); 
            

            ZI = griddata(targetL(:,1),targetL(:,2),targetL(:,3),XI,YI);

%             Tes = delaunay3(targetL(:,1),targetL(:,2),targetL(:,3));
            record=[];
            position=[];
%             k=1;
            sum_CI=zeros(size(XI));
            winPoints=13;
            constant_time=155-1*winPoints;
            constant_time=0;
            for timenow=1:winPoints
                CI = griddata(targetL(:,1),targetL(:,2),targetC(:,timenow+constant_time),XI,YI);
                
%                 C=targetC(:,timenow+constant_time)/max(targetC(:,timenow+constant_time));
%                 CI = griddata(targetL(:,1),targetL(:,2),C,XI,YI);
				TA=isnan(CI);
				ind=find(TA==1);
				CI(ind)=0;

                sum_CI=sum_CI+CI;
                
%                 strength_now=targetC(:,timenow+constant_time);
%                 strength_norm=strength_now/max(strength_now);
%                 ind=find(strength_norm<0.35); %0.5 uAmm as the 
%                 strength_now(ind)=0;
%                 
            end% end of timenow
            if phTime==1
                pre_current=sum_CI./winPoints;
                sum_CI=sum_CI/max(max(sum_CI));
				ind1=find(sum_CI>0.85);
                ActA(k,i)=length(ind1)/size(sum_CI,1)/size(sum_CI,2);
            else
                cur_current=sum_CI./winPoints;
%                 OI(i,phTime-1)=sum(sum(pre_current.*cur_current,1),2)/sum(sum(pre_current,1),2);
                AOI(i,phTime-1)=sum(sum(pre_current.*cur_current,1),2)/size(sum_CI,1)/size(sum_CI,2);
                cur_current=pre_current;
            end
        end %end of phTime
%         current{i}=sum_CI./timePoints;
%         current{i}=sum_CI./5;
    end % end of task
%     C1=current{1};
%     C2=current{2};
%     OI1_2=sum(sum(C1.*C2,1),2)/sum(sum(C2,1),2)
%     OI2_1=sum(sum(C1.*C2,1),2)/sum(sum(C1,1),2)
    
% 	beep
end % end of subject Number

beep
