% This function read the CCD information in the target local area.
function [max_amp,avg_amp]=target_ccd_amp_time_BL(subject)
% function target_ccd (subject,cdr_file,edge_file,timePoints,plotCort)
% function plot_target_ccd (subjectName,cdr_file,edge_file,timePoints,plotCort)
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
winPoints=13;
for k=1:subNo
    subjectName=subject{k}
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
            cdr_file=[cur_task,'_1.cdr'];
            cdr_BL_file=[cur_task,'_BL','.cdr'];
            disp('********************');
            disp(cdr_file)
			cdr_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_file];
			cdr_BL_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_BL_file];

			[cortexL,targetL,targetC]=find_target_PC3(cdr_file_name,edge,timePoints,plotCort);
% 			[cortexL,targetL,targetC]=find_target_PC1(cdr_file_name,edge1_file_name,edge2_file_name,edge3_file_name,edge4_file_name,timePoints,plotCort);
			targetC=targetC(:,end-winPoints+1:end);
            minX=min(targetL(:,1));
            maxX=max(targetL(:,1));
            minY=min(targetL(:,2));
            maxY=max(targetL(:,2));
            minZ=min(targetL(:,3));
            maxZ=max(targetL(:,3));
            x=[minX-5:0.8:maxX+5];
            y=[minY-5:0.8:maxY+5];
            [XI,YI] = meshgrid(x,y); 
            targetC_BL=[];
            [cortexL,targetL,targetC_BL]=find_target_PC3(cdr_BL_file_name,edge,14,plotCort);
            targetC_BL=mean(targetC_BL,2);
            targetC=targetC-repmat(targetC_BL,1,winPoints);


            ZI = griddata(targetL(:,1),targetL(:,2),targetL(:,3),XI,YI);

%             Tes = delaunay3(targetL(:,1),targetL(:,2),targetL(:,3));
%             k=1;
            record=[];
            position=[];
            sum_CI=zeros(size(XI));
            max_C=0;
            constant_time=155-1*winPoints;

            for timenow=1:winPoints
                CI = griddata(targetL(:,1),targetL(:,2),targetC(:,timenow),XI,YI);
				TA=isnan(CI);
				ind=find(TA==1);
				CI(ind)=0;
                [CI_xDim,CI_yDim]=size(CI);

                sum_CI=sum_CI+CI;
                max_C=max_C+max(targetC(:,timenow));
            end% end of timenow
            avg_amp(k,i)=sum(sum(sum_CI,1),2)/(CI_xDim*CI_yDim-length(ind))/winPoints;
            max_amp(k,i)=max_C/winPoints;

        end %end of phTime
    end % end of task
end % end of subject Number
beep
