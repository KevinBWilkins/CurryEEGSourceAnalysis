% This function read the CCD information in the target local area during a time window and then find size of the active area.
% function [ActiveRatio,AOI,LNo,AOI_CI,LNo_CI]=target_ccd_area_BL1(subject)
function [OverArea_r]=target_ccd_OI_time_stat(subject,group)
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


Base_dir=['F:\data\inverse_results\',group,'\'];
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
% LoopNo=fix((timePoints-winOverlap)/(winSize-winOverlap));
LoopNo=1;
needed_len=(LoopNo-1)*(winSize-winOverlap)+winSize;
startPoints=timePoints-needed_len-1;
taskColorMap = ['jet';'jet'];
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
%     winPoints=15;
    
    for i=1:length(taskList)
        
        
        cur_task=taskList{i};
        
        ActCount=0;
%         cdr_file=[cur_task,'2DOF.cdr'];
        cdr_file=[cur_task,'-full.cdr'];
        cdr_BL_file=[cur_task,'_BL','.cdr'];
%             cdr_file=[cur_task,'.cdr'];
        mean_ccd_file=[Base_dir,subjectName,'\new\',cur_task,'.dat'];
        
        disp('********************');
        disp(cdr_file)
        cdr_file_name=[Base_dir,subjectName,'\new\',cdr_file];
        cdr_BL_file_name=[Base_dir,subjectName,'\new\',cdr_BL_file];
% 			[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
        [cortexL,targetL,targetC_wh_time(i,:,:)]=find_target_PC4(cdr_file_name,mask_file_name,timePoints,plotCort);
%         load targetC_wh_time
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
    end %end of i(task) 
    
    for phTime=1:LoopNo+1
        
        if phTime == 1
            cur_window = [1:13];
        else
            cur_window = [timePoints-12:timePoints];
           
        end
        for timeP=1:length(cur_window)
            for i=1:length(taskList)
                cur_task=taskList{i};
                targetC_mean=squeeze(targetC_wh_time(i,:,cur_window(timeP)));
                targetC_mean=targetC_mean';
                [locationNo,TimeNo]=size(targetC_mean);

%                 targetC_mean=targetC;

    %----------save the mean for showing purpose---------------
%             save_data=[save_data,targetC_mean];
%             save (mean_ccd_file,'save_data','-ascii');
    %---------------------------------------------------------------
            

                maxStrenght(k,phTime,i)=max(targetC_mean);
                targetC_norm=targetC_mean/maxStrenght(k,phTime,i);


                [muhat, muci] = expfit(targetC_norm,alpha_static);
                thresh=muci(2);

                ind = find(targetC_norm>thresh);
                active_ind(i,:) = (i+1)*ones (locationNo, 1);
                active_ind(i,ind)= 0;

                ActiveRatio(k,phTime,i)=length(ind)/locationNo; %ActiveRatio: the percentage of actived area

%                 activeC = targetC_mean (ind,:);
                activeL = targetL (ind,:);
%                 CoG(k,phTime,i,:)=sum(repmat(activeC,1,3).*activeL,1)/sum(activeC);
%                 CCD_task(:,i)=targetC_mean;
                CCD_task_norm(:,i)=targetC_norm;

                if debug_yj
                    if cur_task(1)=='a'
                        figure (2) %*phTime)

                    else
                        figure (2*phTime+1)
                    end
                    clf
                    hold on
                    CI = griddata(targetL(:,1),targetL(:,2),targetC_mean,XI,YI);
                    surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
                    h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
                    h1=plot3(CoG(k,phTime,i,1),CoG(k,phTime,i,2),CoG(k,phTime,i,3),'m*','MarkerSize',6);
                    view(90,90)
                    colorbar
                    h=plot3(activeL(:,1),activeL(:,2),activeL(:,3),'mo');
                    title (cur_task);

                end
%                 if phTime~=1 fig_h=plot_data_on_cortex(subjectName, -21, Base_dir, targetL, censul, targetC_norm, cortexL, phTime-1, taskColorMap(i,:), i, cur_window);  end
            end % end of task

            OverArea_norm=CCD_task_norm(:,1).*CCD_task_norm(:,2); %Overlapped area after normalization of the strength
%             OverArea=CCD_task(:,1).*CCD_task(:,2); %Overlapped area without normalization
%             OI(timeP,phTime,1)=sum(OverArea)/sum(CCD_task(:,2)); %Overlap Index: Overlapped area / the total area of the secondary direction
%             OI(timeP,phTime,2)=sum(OverArea)/sum(CCD_task(:,1));

            [muhat, muci] = expfit(OverArea_norm, alpha_static);
            thresh=muci(2);
            ind3=find(OverArea_norm>thresh);
            OverArea_r(timeP+(phTime-1)*13,1)=length(ind3)/locationNo;  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
            
%             OAA(timeP,phTime)=2*length(ind3)/(sum(ActiveRatio(k,phTime,:),3)*locationNo);  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
% 
%             diff=active_ind(1,:)-active_ind(2,:);
%             ind4=find(diff==0);
%             OverArea_r1(timeP,phTime)=length(ind4)/locationNo;  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
%             OAA1(timeP,phTime)=2*length(ind4)/(sum(ActiveRatio(k,phTime,:),3)*locationNo);  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
        end % end of timeP
    end %end of phTime
%     figure
%     boxplot(OverArea_r)
%---------------
%         figure (100)
    win={'1','1','1','1','1','1','1','1','1','1','1','1','1',...
        '2','2','2','2','2','2','2','2','2','2','2','2','2'};
    p=anova1(OverArea_r,win)
    
   

end % end of subject Number



beep
