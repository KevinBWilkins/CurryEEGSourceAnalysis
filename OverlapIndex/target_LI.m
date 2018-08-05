% This function read the CCD information in the target local area and then calculate the laterality index (LI).
% function [ActiveRatio,AOI,LNo,AOI_CI,LNo_CI]=target_ccd_area_BL1(subject)
function [LI,LI_int,LI_act]=target_LI(subject,arm)
% OI is the overlap index taking accout both the area and amplitude
% ActiveRatio is ratio between active location and the total SM area
% OverArea_s is the overlaped active area (area only), no threshould is
%          applied
% AOI_L is the overlaped active location where the activity
%          is larger than 0.25
% AOI_inter is the overlaped area (i.e., obtained based on the interplation
%          on the ective location.) where the activity is larger than 0.25
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
debug_yj=0;    
winPoints=15;
timePoints=155;

% group={};
edgeL1=[];
edgeL2=[];
edgeR1=[];
edgeR2=[];
for k=1:subNo
    subjectName=subject{k};
	edgeL1_file_name=[Base_dir,subjectName,'\Data\MRI\','SM1_L.sp'];
    edgeL2_file_name=[Base_dir,subjectName,'\Data\MRI\','PREMOTOR_L.sp'];
	edgeR1_file_name=[Base_dir,subjectName,'\Data\MRI\','SM1_R.sp'];
    edgeR2_file_name=[Base_dir,subjectName,'\Data\MRI\','PREMOTOR_R.sp'];
	
    [edgeL1,Ecount,ENR]=read_Curry_file3(edgeL1_file_name,'LOCATION',0,0);
	[edgeL2,Ecount,ENR]=read_Curry_file3(edgeL2_file_name,'LOCATION',0,0);
	[edgeR1,Ecount,ENR]=read_Curry_file3(edgeR1_file_name,'LOCATION',0,0);
	[edgeR2,Ecount,ENR]=read_Curry_file3(edgeR2_file_name,'LOCATION',0,0);
	edgeL=[edgeL1;edgeL2];
    edgeR=[edgeR1;edgeR2];
    
    censul_file_name=[Base_dir,subjectName,'\Data\MRI\','CS_L.sp'];
    [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name,'LOCATION',0,0);
    censul_file_name=[Base_dir,subjectName,'\Data\MRI\','CS_R.sp'];
    [censul2,Ecount,ENR]=read_Curry_file3(censul_file_name,'LOCATION',0,0);
    censul=[censul1;censul2];
    record=[];
    for i=1:length(taskList)
        cur_task=taskList{i};
        
        ActCount=0;
        cdr_file=[cur_task,'.cdr'];
        cdr_BL_file=[cur_task,'_BL','.cdr'];
%             cdr_file=[cur_task,'.cdr'];

        disp('********************');
        disp(cdr_file)
		cdr_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_file];
		cdr_BL_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_BL_file];
% 			[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
		[cortexL,targetL_LR{1},targetC_LR{1}]=find_target_PC3(cdr_file_name,edgeL,timePoints,plotCort);
		[cortexL,targetL_LR{2},targetC_LR{2}]=find_target_PC3(cdr_file_name,edgeR,timePoints,plotCort);
        [cortexL,targetL,targetC_BL_LR{1}]=find_target_PC3(cdr_BL_file_name,edgeL,14,plotCort);
        [cortexL,targetL,targetC_BL_LR{2}]=find_target_PC3(cdr_BL_file_name,edgeR,14,plotCort);

        for site=1:2
            targetL=targetL_LR{site};
            targetC=targetC_LR{site};
            targetC_BL=targetC_BL_LR{site};
            targetC=targetC(:,end-winPoints+1:end);
            targetC_BL=mean(targetC_BL,2);
            
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
            
            [locationNo,TimeNo]=size(targetC);
            
            targetC=targetC-repmat(targetC_BL,1,winPoints);
            targetC_sum=sum(targetC,2);
            targetC_norm=targetC_sum/max(targetC_sum);
            
            CI = griddata(targetL(:,1),targetL(:,2),targetC_sum,XI,YI);
            
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
            
            %calculate the active area ratio without interp
            [muhat, muci] = expfit(targetC_norm);
            ind=find(targetC_norm>muci(2));
            ActiveRatio(k,i,site)=length(ind)/locationNo;
            
            %calcultate the activation ratio without interp
            Activation(k,i,site)=sum(targetC_sum);
            
            %calculate the active area ratio with interp
			TA=isnan(CI);
			ind3=find(TA==1);
			CI(ind3)=0;            
			[counts,xx]=imhist(CI);
		    [muhat, muci] = expfit(counts);
		    ind=find(counts>muci(2));
		    thr=xx(ind(end));
		    ind2=find(CI>thr);
		    ActiveA(k,i,site)=length(ind2);
                
%             CCD_task(:,i,site)=targetC_sum;
%             CCD_task_norm(:,i,site)=targetC_norm;
        end %end of site: Left VS. Right
        switch arm{k}
            case {'r','R'}
                LI(k,i)=(ActiveRatio(k,i,1)-ActiveRatio(k,i,2))/(ActiveRatio(k,i,1)+ActiveRatio(k,i,2));
                LI_int(k,i)=(ActiveA(k,i,1)-ActiveA(k,i,2))/(ActiveA(k,i,1)+ActiveA(k,i,2));
                LI_act(k,i)=(Activation(k,i,1)-Activation(k,i,2))/(Activation(k,i,1)+Activation(k,i,2));
            case {'l','L'}
                LI(k,i)=(ActiveRatio(k,i,2)-ActiveRatio(k,i,1))/(ActiveRatio(k,i,1)+ActiveRatio(k,i,2));
                LI_int(k,i)=(ActiveA(k,i,2)-ActiveA(k,i,1))/(ActiveA(k,i,1)+ActiveA(k,i,2));  
                LI_act(k,i)=(Activation(k,i,2)-Activation(k,i,1))/(Activation(k,i,1)+Activation(k,i,2));
        end
    end % end of task
%     OverArea_norm=CCD_task_norm(:,1).*CCD_task_norm(:,2);
%     OverArea=CCD_task(:,1).*CCD_task(:,2);
%     OI(k,1)=sum(OverArea)/sum(CCD_task(:,2));
%     OI(k,2)=sum(OverArea)/sum(CCD_task(:,1));
%     CI_task1 = griddata(targetL(:,1),targetL(:,2),CCD_task_norm(:,1),XI,YI);
%     CI_task2 = griddata(targetL(:,1),targetL(:,2),CCD_task_norm(:,2),XI,YI);
%     TA=isnan(CI_task1);
% 	ind3=find(TA==1);
% 	CI_task1(ind3)=0;
%     TA=isnan(CI_task2);
% 	ind3=find(TA==1);
% 	CI_task2(ind3)=0;
%    
%     OverArea_inter= CI_task1.* CI_task2;
% 
%     OI1(k,1)=sum(sum(OverArea_inter,1),2)/sum(sum(CI_task2,1),2);
%     OI1(k,2)=sum(sum(OverArea_inter,1),2)/sum(sum(CI_task1,1),2);
% 
%     [muhat, muci] = expfit(OverArea_norm);
%     ind1=find(targetC_norm>muci(2));
%     AOI(k)=length(ind1);
%     LNo(k)=locationNo;
%     AOI_L(k)=AOI(k)/LNo(k);
%     OverArea_s(k)=sum(OverArea_norm)/locationNo;    
% 
%     
%     if show_frame
%         figure 
%         hold on
%         surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
%         h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
%         colorbar
%         set(h,'LineWidth',7)
%         view(90,90)
%     end
% 
%     CI = griddata(targetL(:,1),targetL(:,2),OverArea_norm,XI,YI);
% 	TA=isnan(CI);
% 	ind3=find(TA==1);
% 	CI(ind3)=0;
% 
% 	[counts,xx]=imhist(CI);
%     [muhat, muci] = expfit(counts);
%     ind=find(counts>muci(2));
%     thr=xx(ind(end));
%     ind2=find(CI>thr);
%     AOI_CI(k)=length(ind2);
%     ind4=find(CI>0);
%     LNo_CI(k)=length(ind4);
%     AOI_inter(k)=AOI_CI(k)/LNo_CI(k);

%     clear CCD_task CCD_task_norm OverArea OverArea_norm XI YI ZI CI
end % end of subject Number

beep
