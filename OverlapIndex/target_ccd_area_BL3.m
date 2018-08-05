% This function read the CCD information in the target local area and then find size of the active area.
% function [ActiveRatio,AOI,LNo,AOI_CI,LNo_CI]=target_ccd_area_BL1(subject)
function [OI,ActiveRatio,OverArea_s,maxStrenght, CoG]=target_ccd_area_BL3(subject)
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


Base_dir='F:\data\inverse_results\stroke-post\';
subNo=length(subject);
taskList={'abd';'ef'};
phaseList={'whole'};

plotCort=1;

show_frame=0;
yjcolor=['b','r','m'];
currents={};
V1='BL ';
V2='EXE';
debug_yj=1;
% group={};
for k=1:subNo
    subjectName=subject{k};   
    %--------------------------------------------
    %Read in central sulcus file, may be different name
    censul_file_name1=[Base_dir,subjectName,'\Points\','CS_cc.pom'];   %central sulcus
    mask_file_name=[Base_dir,subjectName,'\Points\','Mask.mat'];   %index of ROI
    
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
        
        ActCount=0;
        cdr_file=[cur_task,'.cdr'];
        cdr_BL_file=[cur_task,'_BL','.cdr'];
%             cdr_file=[cur_task,'.cdr'];

        disp('********************');
        disp(cdr_file)
        cdr_file_name=[Base_dir,subjectName,'\new\',cdr_file];
        cdr_BL_file_name=[Base_dir,subjectName,'\new\',cdr_BL_file];
% 			[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
        [cortexL,targetL,targetC]=find_target_PC4(cdr_file_name,mask_file_name,timePoints,plotCort);

        if (i==1 )
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
        [cortexL,targetL,targetC_BL]=find_target_PC4(cdr_BL_file_name,mask_file_name,14,plotCort);
        targetC_BL=mean(targetC_BL,2);
% 			[cortexL,targetC]=find_strength_whole(cdr_file_name,timePoints);
% 			[cortexL,targetC_BL]=find_strength_whole(cdr_BL_file_name,timePoints);

        [locationNo,TimeNo]=size(targetC);

        targetC=targetC-repmat(targetC_BL,1,winPoints);
        targetC_mean=mean(targetC,2);
        maxStrenght(k,i)=max(targetC_mean);
        targetC_norm=targetC_mean/maxStrenght(k,i);
        CoG(k,i,:)=sum(repmat(targetC_mean,1,3).*targetL,1)/sum(targetC_mean);

        if debug_yj
            figure
            hold on
            CI = griddata(targetL(:,1),targetL(:,2),targetC_mean,XI,YI);
            surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
            h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
            h1=plot3(CoG_H(i,1),CoG_H(i,2),CoG_H(i,3),'m*','MarkerSize',6);
            view(90,90)
            colorbar

%                 pause
        end

        ind=find(targetC_norm<0);
%         figure; plot(targetC_norm(ind));
        targetC_norm_abs=abs(targetC_norm);
        [muhat, muci] = expfit(targetC_norm_abs);
        ind=find(targetC_norm_abs>muci(1));

%             ind=find(targetC_norm>0.75);
        ActiveRatio(k,i)=length(ind)/locationNo; %ActiveRatio: the percentage of actived area
        CCD_task(:,i)=targetC_mean;
        CCD_task_norm(:,i)=targetC_norm;
    end % end of task
    OverArea_norm=abs(CCD_task_norm(:,1).*CCD_task_norm(:,2)); %Overlapped area after normalization of the strength
    OverArea=abs(CCD_task(:,1).*CCD_task(:,2)); %Overlapped area without normalization
    OI(k,1)=sum(OverArea)/sum(abs(CCD_task(:,2))); %Overlap Index: Overlapped area / the total area of the secondary direction
    OI(k,2)=sum(OverArea)/sum(abs(CCD_task(:,1)));
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
%     OI_inter(k,1)=sum(sum(OverArea_inter,1),2)/sum(sum(CI_task2,1),2); %Overlap index: Overlapped area / the total area of the secondary direction; Differently, first interplate data to a grid
%     OI_inter(k,2)=sum(sum(OverArea_inter,1),2)/sum(sum(CI_task1,1),2);

    [muhat, muci] = expfit(OverArea_norm);
    ind1=find(targetC_norm>muci(1));
    OverArea_s(k)=length(ind1)/locationNo;  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area

    
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

    clear CCD_task CCD_task_norm OverArea OverArea_norm XI YI ZI CI
end % end of subject Number

beep
