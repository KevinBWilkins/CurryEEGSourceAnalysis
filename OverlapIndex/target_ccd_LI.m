% This function read the CCD information in the target local area and then find size of the active area.
% function [ActiveRatio,AOI,LNo,AOI_CI,LNo_CI]=target_ccd_area_BL1(subject)
%function [OI,ActiveRatio,OverArea_r,OAA, OverArea_r1,OAA1, maxStrenght, CoG_s]=target_ccd_LI(subject,group,startTime)
function [LI,targetC_wh_time_ROI]=target_ccd_LI(subject)
% OI is the overlap index taking accout both the area and amplitude
% ActiveRatio is ratio between active location and the total SM area
% OverArea_s is the overlaped active area (area only), no threshould is
%          applied
% AOI_L is the overlaped active location where the activity
%          is larger than 0.25
% AOI_inter is the overlaped area (i.e., obtained based on the interplation
%          on the ective location.) where the activity is larger than 0.25
%method={'1'};
% pre_method=method(1);
%addpath('E:\Matlab\Eeg\inv')


%Base_dir=['F:\data\inverse_results\',group,'\'];
%Base_dir='C:\Users\Kevin\Desktop\ImagingR01\';
%Base_dir='C:\Users\Kevin\Desktop\LP=1_Results\';
%Base_dir='D:\LP=1_Results\'
Base_dir='D:\sLORETA\Stroke\'
subNo=length(subject);
taskList={'sLORETA_LiftOpen'};


plotCort=1;

show_frame=1;
debug_yj=0;

for k=1:subNo
    subjectName=subject{k};   
    %--------------------------------------------
    %Read in central sulcus file, may be different name
    censul_file_name1=[Base_dir,subjectName,'\Points\','CS-KW.pom'];   %central sulcus
    mask_file_name=[Base_dir,subjectName,'\Points\','Mask.mat'];   %index of ROI
    
    %if only one central sulcus file
    [censul,Ecount,ENR]=read_Curry_file3(censul_file_name1,'LOCATION',0,0);

   
    %----------------------------------------------------------------------

    for i=1:length(taskList)
        
        cur_task=taskList{i};
        timePoints=14;
        

        cdr_file=[cur_task,'.cdr'];
        %cdr_BL_file=['lift day2.cdr'];
%             cdr_file=[cur_task,'.cdr'];
       % mean_ccd_file=[Base_dir,subjectName,'\new\',cur_task,'.dat']; %????
        
        disp('********************');
        disp(cdr_file)
        cdr_file_name=[Base_dir,subjectName,'\new\',cdr_file];
        %cdr_BL_file_name=[Base_dir,subjectName,'\new\',cdr_BL_file];
% 			[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
        %[cortexL,targetL,targetC]=find_target_PC4(cdr_file_name,mask_file_name,timePoints,plotCort);
        [cortexL,targetL_cell,targetC_wh_time_ROI(:,i)]=find_target_PC4(cdr_file_name,mask_file_name,timePoints,plotCort,0:1);

        if (i==1 )
            targetL = cell2mat(targetL_cell);
            
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
        
       
        target_right{i}=targetC_wh_time_ROI{1,i};
        mTarget_right{i}=mean(target_right{i},2);
        target_left{i}=targetC_wh_time_ROI{2,i};
        mTarget_left{i}=mean(target_left{i},2);
        target_all_activity = [mTarget_right{i};mTarget_left{i}];
        threshold = .25* max (target_all_activity);  % Set threshold. Set to .01 for no threshold
        tmp= mTarget_left{i};
        ind_less_threshold = find (tmp <threshold);
        tmp(ind_less_threshold) = 0;
        mTarget_left{i} = tmp;
        tmp= mTarget_right{i};
        ind_less_threshold = find (tmp <threshold);
        tmp(ind_less_threshold) =0;
        mTarget_right{i} = tmp;
        
        
        Contral = mean (mTarget_left{i});  % Contral=contralateral hemisphere; right hand affected=make this left hemisphere
        Ipsi = mean (mTarget_right{i});
        LI(k,i) = (Contral-Ipsi)./(Contral+Ipsi);
        
        
        if show_frame
            figure 
            target_all_activity = [mTarget_right{i};mTarget_left{i}];
            CI = griddata(targetL(:,1),targetL(:,2),target_all_activity,XI,YI);
            hold on
            surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
            h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
            colorbar
            set(h,'LineWidth',7)
            view(90,90)
        end

        
        
    end %end of i_task
    
%---------------------------------------------------------------
    
    if debug_yj
%         overlappedL = targetL (ind3,:);
%         overlappedL1 = targetL (ind4,:);
%         subject_pic_Dir=[Base_dir,subjectName,'\pic'];
%         if ~exist(subject_pic_Dir,'dir')
%             mkdir(subject_pic_Dir);
%         end
%         for Pic_N=1:2;
%             if Pic_N==1 
%                 cur_task=taskList{1};
%                 figure (2)
%                 h=plot3(overlappedL(:,1),overlappedL(:,2),overlappedL(:,3),'m.');
%             else
%                 cur_task=taskList{2};
%                 figure (4)
%                 h=plot3(overlappedL1(:,1),overlappedL1(:,2),overlappedL1(:,3),'m.');
%             end
%             
%             fig_file_name=[Base_dir,subjectName,'\pic\',subjectName,cur_task,'.jpg'];   %index of ROI
%             saveas (h, fig_file_name, 'jpg');
%         end

    end

    
end % end of subject Number

beep
