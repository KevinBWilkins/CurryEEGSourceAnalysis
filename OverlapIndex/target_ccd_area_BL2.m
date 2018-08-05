% This function read the CCD information in the target local area and then find size of the active area.
% function [ActiveRatio,AOI,LNo,AOI_CI,LNo_CI]=target_ccd_area_BL1(subject)
function [OI,OI1,ActiveRatio,OverArea_s,AOI_L,AOI_inter]=target_ccd_area_BL2(subject)
% OI is the overlap index taking accout both the area and amplitude
% ActiveRatio is ratio between active location and the total SM area
% OverArea_s is the overlaped active area (area only), no threshould is
%          applied
% AOI_L is the overlaped active location where the activity
%          is larger than 0.25
% AOI_inter is the overlaped area (i.e., obtained based on the interplation
%          on the ective location.) where the activity is larger than 0.25
method={'1'};
% pre_method=method(1);
addpath('C:\Documents and Settings\Jun Yao\Matlab\Eeg\inv')


Base_dir='F:\data\inverse_results\';
subNo=length(subject);
taskList={'abd';'ef'};
phaseList={'whole'};

plotCort=0;

show_frame=0;
yjcolor=['b','r','m'];
currents={};
V1='BL ';
V2='EXE';
debug_yj=1;
% group={};
for k=1:subNo
    subjectName=subject{k};
    edge1_file_name=[Base_dir,subjectName,'\Points\','M1_rt_cc.pom'];%M1_rt
    edge2_file_name=[Base_dir,subjectName,'\Points\','M1_lt_cc.pom'];%M1_lt
    edge3_file_name=[Base_dir,subjectName,'\Points\','S1_rt_cc.pom'];%S1_rt
    edge4_file_name=[Base_dir,subjectName,'\Points\','S1_lt_cc.pom'];%S1_lt
	edge5_file_name=[Base_dir,subjectName,'\Points\','PM_rt_cc.pom'];%PM_rt
    edge6_file_name=[Base_dir,subjectName,'\Points\','PM_lt_cc.pom'];%PM_lt
    edge7_file_name=[Base_dir,subjectName,'\Points\','SMA_rt_cc.pom'];%SMA_rt
    edge8_file_name=[Base_dir,subjectName,'\Points\','SMA_lt_cc.pom'];%SMA_lt
	
    [edge1,Ecount,ENR]=read_Curry_file4_AC(edge1_file_name,'LOCATION',0,0);
	[edge2,Ecount,ENR]=read_Curry_file4_AC(edge2_file_name,'LOCATION',0,0);
	[edge3,Ecount,ENR]=read_Curry_file4_AC(edge3_file_name,'LOCATION',0,0);
	[edge4,Ecount,ENR]=read_Curry_file4_AC(edge4_file_name,'LOCATION',0,0);
    [edge5,Ecount,ENR]=read_Curry_file4_AC(edge5_file_name,'LOCATION',0,0);
	[edge6,Ecount,ENR]=read_Curry_file4_AC(edge6_file_name,'LOCATION',0,0);
	[edge7,Ecount,ENR]=read_Curry_file4_AC(edge7_file_name,'LOCATION',0,0);
	[edge8,Ecount,ENR]=read_Curry_file4_AC(edge8_file_name,'LOCATION',0,0);

    edge = [edge1;edge2;edge3;edge4;edge5;edge6;edge7;edge8];
    areanum = 1;
% 	edge=[edge1;edge3];
    
    %--------------------------------------------
    %Read in central sulcus file, may be different name
    censul_file_name1=[Base_dir,subjectName,'\Points\','CS_cc.pom'];   %central sulcus
    
    %if only one central sulcus file
    [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name1,'LOCATION',0,0);
    censul2 = censul1;
    censul=[censul1];
    %----------------------------------------------------------------------

    record=[];
    winPoints=15;
    for i=1:length(taskList)
        cur_task=taskList{i};
        timePoints=155;
        
        for methodInd=1:length(method)
            ActCount=0;
            cdr_file=[cur_task,'.cdr'];
            cdr_BL_file=[cur_task,'_BL','.cdr'];
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
                CI = griddata(targetL(:,1),targetL(:,2),targetC_sum/winPoints,XI,YI);
                surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
                h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
                view(90,90)
                colorbar
                
%                 pause
            end
            
            [muhat, muci] = expfit(targetC_norm);
            ind=find(targetC_norm>muci(2));
            
%             ind=find(targetC_norm>0.75);
            ActiveRatio(k,i)=length(ind)/locationNo; %ActiveRatio: the percentage of actived area
            
        end %end of methodInd
        CCD_task(:,i)=targetC_sum;
        CCD_task_norm(:,i)=targetC_norm;
    end % end of task
    OverArea_norm=CCD_task_norm(:,1).*CCD_task_norm(:,2); %Overlapped area after normalization of the strength
    OverArea=CCD_task(:,1).*CCD_task(:,2); %Overlapped area without normalization
    OI(k,1)=sum(OverArea)/sum(CCD_task(:,2)); %Overlap Index: Overlapped area / the total area of the secondary direction
    OI(k,2)=sum(OverArea)/sum(CCD_task(:,1));
    CI_task1 = griddata(targetL(:,1),targetL(:,2),CCD_task_norm(:,1),XI,YI);
    CI_task2 = griddata(targetL(:,1),targetL(:,2),CCD_task_norm(:,2),XI,YI);
    TA=isnan(CI_task1);
	ind3=find(TA==1);
	CI_task1(ind3)=0;
    TA=isnan(CI_task2);
	ind3=find(TA==1);
	CI_task2(ind3)=0;
   
    OverArea_inter= CI_task1.* CI_task2;

    OI1(k,1)=sum(sum(OverArea_inter,1),2)/sum(sum(CI_task2,1),2); %Overlap index: Overlapped area / the total area of the secondary direction; Differently, first interplate data to a grid
    OI1(k,2)=sum(sum(OverArea_inter,1),2)/sum(sum(CI_task1,1),2);

    [muhat, muci] = expfit(OverArea_norm);
    ind1=find(targetC_norm>muci(2));
    AOI(k)=length(ind1);
    LNo(k)=locationNo;
    AOI_L(k)=AOI(k)/LNo(k);
    OverArea_s(k)=sum(OverArea_norm)/locationNo;  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area

    
    if show_frame
        figure 
        CI = griddata(targetL(:,1),targetL(:,2),OverArea,XI,YI);
        hold on
        surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
        h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
        colorbar
        set(h,'LineWidth',7)
        view(90,90)
    end

    CI = griddata(targetL(:,1),targetL(:,2),OverArea_norm,XI,YI);
	TA=isnan(CI);
	ind3=find(TA==1);
	CI(ind3)=0;

	[counts,xx]=imhist(CI);
    [muhat, muci] = expfit(counts);
    ind=find(counts>muci(2));
    thr=xx(ind(end));
    ind2=find(CI>thr);
    AOI_CI(k)=length(ind2);
    ind4=find(CI>0);
    LNo_CI(k)=length(ind4);
    AOI_inter(k)=AOI_CI(k)/LNo_CI(k); % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area with the data interplate to the grid first

    clear CCD_task CCD_task_norm OverArea OverArea_norm XI YI ZI CI
end % end of subject Number

beep
