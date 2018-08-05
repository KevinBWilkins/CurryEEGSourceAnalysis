% This function read the CCD information in the target local area and then find size of the active area.
% function [ActiveRatio,AOI,LNo,AOI_CI,LNo_CI]=target_ccd_area_BL1(subject)
function [ActiveRatio,OverArea_s,AOI_L,AOI_inter]=target_ccd_area_BL1(subject)

method={'1'};
% pre_method=method(1);
addpath('C:\Documents and Settings\Jun Yao\Matlab\Eeg\inv')


Base_dir='C:\Documents and Settings\Jun Yao\Subject\';
subNo=length(subject);
taskList={'abd';'flex'};
phaseList={'whole'};

plotCort=0;

show_frame=1;
yjcolor=['b','r','m'];
currents={};
V1='BL ';
V2='EXE';
debug_yj=0;
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
% 	edge=[edge1;edge3];
    
    censul_file_name=[Base_dir,subjectName,'\Data\MRI\','CS_L.sp'];
    [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name,'LOCATION',0,0);
    censul_file_name=[Base_dir,subjectName,'\Data\MRI\','CS_R.sp'];
    [censul2,Ecount,ENR]=read_Curry_file3(censul_file_name,'LOCATION',0,0);
    censul=[censul1;censul2];
    record=[];
    winPoints=13;
    for i=1:length(taskList)
        cur_task=taskList{i};
        timePoints=13;
        
        for methodInd=1:length(method)
            ActCount=0;
            cdr_file=[cur_task,'.cdr'];
            cdr_BL_file=[cur_task,'_BL','.cdr'];
            cdr_file=[cur_task,'_',method{methodInd},'.cdr'];            
	%             cdr_file=[cur_task,'.cdr'];
	
            disp('********************');
            disp(cdr_file)
			cdr_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_file];
			cdr_BL_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_BL_file];
	% 			[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
			[cortexL,targetL,targetC]=find_target_PC3(cdr_file_name,edge,timePoints,plotCort);
            
            if (i==1 & methodInd==1)
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
            end
            targetC=targetC(:,end-winPoints+1:end);
            targetC_BL=[];
            [cortexL,targetL,targetC_BL]=find_target_PC3(cdr_BL_file_name,edge,14,plotCort);
            targetC_BL=mean(targetC_BL,2);
% 			[cortexL,targetC]=find_strength_whole(cdr_file_name,timePoints);
% 			[cortexL,targetC_BL]=find_strength_whole(cdr_BL_file_name,timePoints);
            
            [locationNo,TimeNo]=size(targetC);
            
            targetC=targetC-repmat(targetC_BL,1,winPoints);
            targetC_sum=sum(targetC,2);
            targetC_norm=targetC_sum/max(targetC_sum);
            
            
            if debug_yj
                figure
                hold on
                CI = griddata(targetL(:,1),targetL(:,2),targetC_norm,XI,YI);
                surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
                h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
                view(90,90)
                colorbar
                
                pause
            end
            
            ind=find(targetC_norm>0.75);
            ActiveRatio(k,i)=length(ind)/locationNo;
            
        end %end of methodInd
        CCD_task(:,i)=targetC_norm;
    end % end of task

    OverArea=CCD_task(:,1).*CCD_task(:,2);
    ind1=find(OverArea>0.25);
    AOI(k)=length(ind1);
    LNo(k)=locationNo;
    AOI_L(k)=AOI(k)/LNo(k);

    CI = griddata(targetL(:,1),targetL(:,2),OverArea,XI,YI);
    OverArea_s(k)=sum(OverArea)/locationNo;
    
    if show_frame
        figure 
        hold on
        surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
        h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
        colorbar
        set(h,'LineWidth',7)
        view(90,90)
    end
	TA=isnan(CI);
	ind3=find(TA==1);
	CI(ind3)=-1;

    ind2=find(CI>0.25);
    AOI_CI(k)=length(ind2);
    ind4=find(CI>0);
    LNo_CI(k)=length(ind4);
    AOI_inter(k)=AOI_CI(k)/LNo_CI(k);

    clear CCD_task OverArea XI YI ZI CI
end % end of subject Number

beep
