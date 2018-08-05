% This function read the CCD information in the target local area and then find size of the active area.
function ActiveRatio=target_ccd_area_BL(subject)
method={'1'};
% pre_method=method(1);
addpath('C:\Documents and Settings\Jun Yao\Matlab\Eeg\inv')


Base_dir='C:\Documents and Settings\Jun Yao\Subject\';
subNo=length(subject);
taskList={'abd';'flex'};
phaseList={'whole'};

plotCort=0;

show_frame=0;
yjcolor=['b','r','m'];
currents={};
V1='BL ';
V2='EXE';
% group={};
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
        timePoints=13;
        
        for methodInd=1:length(method)
            ActCount=0;
            cdr_file=[cur_task,'_',method{methodInd},'.cdr'];
            cdr_BL_file=[cur_task,'_BL','.cdr'];
	%             cdr_file=[cur_task,'.cdr'];
	
            disp('********************');
            disp(cdr_file)
			cdr_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_file];
			cdr_BL_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_BL_file];
	% 			[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
			[cortexL,targetL,targetC]=find_target_PC3(cdr_file_name,edge,timePoints,plotCort);
            [cortexL,targetL,targetC_BL]=find_target_PC3(cdr_BL_file_name,edge,20,plotCort);
            
% 			[cortexL,targetC]=find_strength_whole(cdr_file_name,timePoints);
% 			[cortexL,targetC_BL]=find_strength_whole(cdr_BL_file_name,timePoints);
            
            [locationNo,TimeNo]=size(targetC);
            for locationInd=1:locationNo
                X=[targetC_BL(locationInd,:),targetC(locationInd,:)];
                group1=repmat(V1,20,1);
                group2=repmat(V2,TimeNo,1);
                cellGroup=cellstr([group1;group2]);
                cellGroup=cellGroup';
                p = anova1(X,cellGroup,'off');
                if p<0.005
                    ActCount=ActCount+1;
                end    
            end %end of locationInd
            ActiveRatio(k,i)=ActCount/locationNo;
        end %end of methodInd
    end % end of task
end % end of subject Number

beep
