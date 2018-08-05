% This function read the CCD information in the target local area during a time window and then find size of the active area.
% function [ActiveRatio,AOI,LNo,AOI_CI,LNo_CI]=target_ccd_area_BL1(subject)
function [act_ratio_loc_time]=target_ccd_OI_timebins_ROI_subLocation(subject,group,arm)
% OI is the overlap index taking accout both the area and amplitude
% ActiveRatio is ratio between active location and the total SM area
% OverArea_s is the overlaped active area (area only), no threshould is
%          applied
% AOI_L is the overlaped active location where the activity
%          is larger than 0.25
% AOI_inter is the overlaped area (i.e., obtained based on the interplation
%          on the ective location.) where the activity is larger than 0.25
%example: [OI,ActiveRatio,OverArea_r,OAA, OverArea_r1,OAA1, maxStrenght,CoG_s]=target_ccd_OI_timebins({'MG'},'stroke-pre');

method={'1'};
% pre_method=method(1);
addpath('E:\Matlab\Eeg\inv')



subNo=length(subject);
taskList={'abd';'ef'};
phaseList={'whole'};

plotCort=1;
plotCoG=1;

show_frame=0;
yjcolor=['b','r','m'];
currents={};
V1='BL ';
V2='EXE';
debug_yj=0;
% group={};
save_data=[];
alpha_static = 0.01;


winOverlap=1;
winSize=13;
timePoints=155;
LoopNo=fix((timePoints-winOverlap)/(winSize-winOverlap));
needed_len=(LoopNo-1)*(winSize-winOverlap)+winSize;
startPoints=timePoints-needed_len-1;
taskColorMap = ['jet';'jet'];
for k=1:subNo
    subjectArm=arm{k}; 
    if subjectArm == 'l'
        location_order = [3 1 7 5 6 8 2 4];
    elseif subjectArm == 'r'
        location_order = [6 8 2 4 3 1 7 5];   
    end
    
    subjectName=subject{k};   
    
    for i_group = 1:2
        cur_group = group{i_group};
        Base_dir=['F:\data\inverse_results\',cur_group,'\'];
        %--------------------------------------------
        %Read in central sulcus file, may be different name
        censul_file_name1=[Base_dir,subjectName,'\Points\','CS_cc.pom'];   %central sulcus
         %mask_file_name=[Base_dir,subjectName,'\Points\','Mask.mat'];   %index of ROI
         mask_file_name=['F:\data\inverse_results\stroke-pre\',subjectName,'\Points\','Mask.mat'];   %index of ROI



        %if only one central sulcus file
        [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name1,'LOCATION',0,0);
        censul2 = censul1;
        censul=[censul1];
        %----------------------------------------------------------------------

        record=[];

        for i=1:length(taskList)


            cur_task=taskList{i};

            ActCount=0;
            % cdr_file=[cur_task,'2DOF.cdr'];
            cdr_file=[cur_task,'-full.cdr'];
            cdr_BL_file=[cur_task,'_BL','.cdr'];
            %  cdr_file=[cur_task,'.cdr'];
            mean_ccd_file=[Base_dir,subjectName,'\new\',cur_task,'.dat'];

            disp('********************');
            disp(cdr_file)
            cdr_file_name=[Base_dir,subjectName,'\new\',cdr_file]
            cdr_BL_file_name=[Base_dir,subjectName,'\new\',cdr_BL_file];
            %	[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
            [cortexL,targetL_all,targetC_wh_time_ROI(:,i)]=find_target_PC4(cdr_file_name,mask_file_name,timePoints,plotCort,1:8);
            targetC_all = cell2mat(squeeze(targetC_wh_time_ROI(:,i)));
            max_targetC(i) = max(max(targetC_all));
            all_targetC(i) = sum (sum(targetC_all));
            total_location_no = size(targetC_all,1);
            % load targetC_wh_time
            if (i==1 )
                    all_targetL = cell2mat(targetL_all);

                    minX=min(all_targetL(:,1));
                    maxX=max(all_targetL(:,1));
                    minY=min(all_targetL(:,2));
                    maxY=max(all_targetL(:,2));
                    minZ=min(all_targetL(:,3));
                    maxZ=max(all_targetL(:,3));
                    x=[minX-5:0.8:maxX+5];
                    y=[minY-5:0.8:maxY+5];
                    [XI,YI] = meshgrid(x,y); 
                    ZI = griddata(all_targetL(:,1),all_targetL(:,2),all_targetL(:,3),XI,YI);
            end
        end %end of i(task) 

        for phTime=1:LoopNo+1

            if phTime <= LoopNo
                cur_window = (phTime-1)*(winSize-winOverlap)+1+startPoints:(phTime-1)*(winSize-winOverlap)+winSize+startPoints;
            else
                cur_window = (phTime-1)*(winSize-winOverlap)+1+startPoints:timePoints;

            end
            for i_location=1:8 % location
                for i=1:length(taskList)
                    cur_task=taskList{i};
                    targetC_cell=targetC_wh_time_ROI{location_order(i_location),i};
                    targetC=targetC_cell(:,cur_window);
                    targetL=targetL_all{location_order(i_location)};

                    [locationNo,TimeNo]=size(targetC);


                    act_ratio_loc_time (k,i,i_location,phTime)= sum(sum(targetC))/all_targetC(i);

                end % end of task
            end %end of the location

        end %end of time bins
        clear XI YI ZI CI

        for i=1:length(taskList)
            figure
            surf(squeeze(act_ratio_loc_time(k,i,:,:)));
        end
    

    
    end % end of group

end % end of subject Number

for i=1: 2 %two tasks
    task_act_loc_time=squeeze(act_ratio_loc_time(:,i,:,:));
    mean_task_loc_time=mean(task_act_loc_time,1);
    mean_task_loc_time=squeeze(mean_task_loc_time); 
    figure
    hold on
    
    if i==1
        title ('SABD');
    else
        title ('EF');
    end
    ste_task_loc_time=std(task_act_loc_time,0,1)./sqrt(7-1);
    ste_task_loc_time=squeeze(ste_task_loc_time);
    surf(ste_task_loc_time+mean_task_loc_time,'FaceAlpha',0.2)
    shading interp
    surfc(mean_task_loc_time)
    set(gca,'YTickLabel',{'cS1','cM1','cPM','cSMA','iSMA','iPM','iM1','iS1'})
end
save (strcat('E:\Documents\paper\overlap-time\',group,'-m' ), 'act_ratio_loc_time')
beep
