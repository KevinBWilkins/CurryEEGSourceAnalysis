% This function read the CCD information in the target local area during a time window and then find size of the active area.
% function [ActiveRatio,AOI,LNo,AOI_CI,LNo_CI]=target_ccd_area_BL1(subject)
function [ActiveRatio, maxStrenght]=target_ccd_OI_timebins_ROI_KW(subject,group,arm)
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
%addpath('E:\Matlab\Eeg\inv')


Base_dir='C:\Users\Kevin\Desktop\LP=1_Results\';
subNo=length(subject);
taskList={'Post-Open-Long'};
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
    subjectName=subject{k};   
    if arm{k} == 'l'
        location_order = [3 1 7 5 6 8 2 4];
    elseif arm{k} == 'r'
        location_order = [6 8 2 4 3 1 7 5];   
    end

    %--------------------------------------------
    %Read in central sulcus file, may be different name
    censul_file_name1=[Base_dir,subjectName,'\Points\','CS-KW.pom'];   %central sulcus
    mask_file_name=[Base_dir,subjectName,'\Points\','Mask.mat'];   %index of ROI
%     mask_file_name=['F:\data\inverse_results\stroke-pre\',subjectName,'\Points\','Mask.mat'];   %index of ROI
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
%         cdr_file=[cur_task,'-full.cdr'];
        cdr_BL_file=[cur_task,'_BL','.cdr'];
            cdr_file=[cur_task,'.cdr'];
        mean_ccd_file=[Base_dir,subjectName,'\new\',cur_task,'.dat'];
        
        disp('********************');
        disp(cdr_file)
        cdr_file_name=[Base_dir,subjectName,'\new\',cdr_file];
        cdr_BL_file_name=[Base_dir,subjectName,'\new\',cdr_BL_file];
% 			[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
        [cortexL,targetL_cell,targetC_wh_time_ROI(:,i)]=find_target_PC4(cdr_file_name,mask_file_name,timePoints,plotCort,1:8);
        all_targetC (i, :,:) = cell2mat(squeeze(targetC_wh_time_ROI(:,i)));
        max_targetC(i) = max(max(squeeze(all_targetC(i,:,:))));
%         load targetC_wh_time
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
    end %end of i(task) 
    
    for phTime=1:LoopNo+1
        
        if phTime <= LoopNo
            cur_window = (phTime-1)*(winSize-winOverlap)+1+startPoints:(phTime-1)*(winSize-winOverlap)+winSize+startPoints;
        else
            cur_window = (phTime-1)*(winSize-winOverlap)+1+startPoints:timePoints;
           
        end
%         for i_location=1:8 % location
            for i=1:length(taskList)
                cur_task=taskList{i};
%                 targetC_cell=targetC_wh_time_ROI{location_order(i_location),i};
%                 targetC=targetC_cell(:,cur_window);
                targetC = squeeze (all_targetC (i, :,cur_window));

                [locationNo,TimeNo]=size(targetC);

                targetC_mean=mean(targetC,2);

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

                activeC = targetC_mean (ind,:);
                activeL = targetL (ind,:);
                CoG(k,phTime,i,:)=sum(repmat(activeC,1,3).*activeL,1)/sum(activeC);
                CCD_task(:,i)=targetC_mean;
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
                if phTime~=1 fig_h=plot_data_on_cortex(subjectName, -21, Base_dir, targetL, censul, targetC_norm, cortexL, phTime-1, taskColorMap(i,:), i, cur_window);  end
            end % end of task
%         end %end of the location
%         OverArea_norm=CCD_task_norm(:,1).*CCD_task_norm(:,2); %Overlapped area after normalization of the strength
%         OverArea=CCD_task(:,1).*CCD_task(:,2); %Overlapped area without normalization
%         OI(k,phTime,1)=sum(OverArea)/sum(CCD_task(:,2)); %Overlap Index: Overlapped area / the total area of the secondary direction
%         OI(k,phTime,2)=sum(OverArea)/sum(CCD_task(:,1));
% 
%         [muhat, muci] = expfit(OverArea_norm, alpha_static);
%         thresh=muci(2);
%         ind3=find(OverArea_norm>thresh);
%         OverArea_r(k,phTime)=length(ind3)/locationNo;  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
%         OAA(k,phTime)=2*length(ind3)/(sum(ActiveRatio(k,phTime,:),3)*locationNo);  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
% 
%         diff=active_ind(1,:)-active_ind(2,:);
%         ind4=find(diff==0);
%         OverArea_r1(k,phTime)=length(ind4)/locationNo;  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
%         OAA1(k,phTime)=2*length(ind4)/(sum(ActiveRatio(k,phTime,:),3)*locationNo);  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
%---------------
%         figure (100)
%         
%         subplot(4,6,1+2*(phTime-1));
%         plot(OverArea_norm); hold on
%         plot(thresh*ones(size(OverArea_norm)),'r')
        

        if debug_yj
            overlappedL = targetL (ind3,:);
            overlappedL1 = targetL (ind4,:);
            subject_pic_Dir=[Base_dir,subjectName,'\pic'];
            if ~exist(subject_pic_Dir,'dir')
                mkdir(subject_pic_Dir);
            end
            for Pic_N=1:2;
                if Pic_N==1 
                    cur_task=taskList{1};
                    figure (2)
                    h=plot3(overlappedL(:,1),overlappedL(:,2),overlappedL(:,3),'m.');
                else
                    cur_task=taskList{2};
                    figure (2*phTime+1)
                    h=plot3(overlappedL1(:,1),overlappedL1(:,2),overlappedL1(:,3),'m.');
                end

                fig_file_name=[Base_dir,subjectName,'\pic\',subjectName,cur_task,'_try.jpg'];   %index of ROI
                saveas (h, fig_file_name, 'jpg');
            end
        end
%         overlappedL = targetL (ind3,:);
%         overlappedC = targetC (ind3,1);
%         plot_data_on_cortex(subjectName, -21, Base_dir, overlappedL, censul, overlappedC, cortexL, phTime, taskColorMap(i,:), 3, cur_window);
% 
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
    %     CoG_s=['(',num2str(CoG(:,1,1)), ', ',num2str(CoG(:,1,2)),', ',num2str(CoG(:,1,3)),')';'(',num2str(CoG(:,2,1)),', ',num2str(CoG(:,2,2)),', ',num2str(CoG(:,2,3)),')'];
%        CoG_s{phTime}={num2str(CoG(:,phTime,1,1)),num2str(CoG(:,phTime,1,2)),num2str(CoG(:,phTime,1,3)),num2str(CoG(:,phTime,2,1)),num2str(CoG(:,phTime,2,2)),num2str(CoG(:,phTime,2,3))};
        clear CCD_task CCD_task_norm OverArea OverArea_norm 
    end %end of time bins
    clear XI YI ZI CI
    
    
    fig_file_name=[Base_dir,subjectName,'\pic\',subjectName,'cortex+final_try.jpg'];   %index of ROI
    saveas (fig_h, fig_file_name, 'jpg');
    
    
    
    %----------calculate distance from CoG(end) to CoG(baseline)------------
    NoBaseline=5;
    
    for i_task=1:length(taskList)
        CoG_end = squeeze (CoG(:,LoopNo,i_task,:));
        i_task
        CoG_dis_tmp=[];
        for i_search=1:NoBaseline
            cur_CoG=squeeze (CoG(:,i_search,i_task,:));
            CoG_dis_tmp(i_search)=sqrt(sum ((CoG_end-cur_CoG).^2));
        end
        
        CoG_dis(i_task) = min(CoG_dis_tmp);
        %----------calculate the mean CoG during three windows [1:4],[5:8];[9:12]CoG(baseline)------------
        A=mean ( squeeze(CoG(:,[1:4],i_task,:)), 1);
        B=mean ( squeeze(CoG(:,[5:8],i_task,:)), 1);
        C=mean ( squeeze(CoG(:,[9:12],i_task,:)), 1);
        mean_CoG = [ A,B,C]
    end
    %------------------
    
    if plotCoG
        for i=1:length(taskList)
            cur_task=taskList{i};
            figure (50+i)
            p=plot3(targetL(:,1),targetL(:,2),targetL(:,3),'LineStyle','none','Marker','.','MarkerFaceColor','b','EraseMode','none', 'MarkerSize',1);
            
            Tmax=(ceil(100*1.1*max([min(targetL(:,1)) max(targetL(:,1)) min(targetL(:,2)) max(targetL(:,2)) min(targetL(:,3)) max(targetL(:,3))]))+0.5)/100;
            Tmin=(ceil(100*1.1*min([min(targetL(:,1)) max(targetL(:,1)) min(targetL(:,2)) max(targetL(:,2)) min(targetL(:,3)) max(targetL(:,3))]))-0.5)/100;
            axis([Tmin Tmax Tmin Tmax Tmin Tmax])
            view ([0 90])
            hold on
            plot3(targetL(:,1),targetL(:,2),targetL(:,3),'LineStyle','none','Marker','.','MarkerFaceColor','b','EraseMode','none', 'MarkerSize',5);
            plot3(censul(:,1),censul(:,2),censul(:,3),'LineStyle','none','Marker','.','MarkerFaceColor','m','EraseMode','none', 'MarkerSize',10);
            for phTime=1:LoopNo
%                 set(p,'XData',CoG(:,phTime,i,1),'YData',CoG(:,phTime,i,2),'ZData',CoG(:,phTime,i,3),'color',[(LoopNo-phTime)/LoopNo phTime/LoopNo phTime/LoopNo],'MarkerSize',5,'LineStyle','none','EraseMode','none')
                plot3(CoG(:,phTime,i,1),CoG(:,phTime,i,2),CoG(:,phTime,i,3),'LineStyle','none','Marker','.','MarkerFaceColor','b','EraseMode','none', 'MarkerSize',20,'color',[(LoopNo-phTime)/LoopNo phTime/LoopNo phTime/LoopNo]);
                drawnow
                hold on
            end
            view ([0 90])
            I=getframe(gcf);
            subject_pic_Dir=[Base_dir,subjectName,'\pic'];
            if ~exist(subject_pic_Dir,'dir')
                mkdir(subject_pic_Dir);
            end

            fig_file_name=[Base_dir,subjectName,'\pic\',subjectName,cur_task,'CoG_try.jpg'];   %index of ROI
            saveas (p, fig_file_name, 'jpg');
        end
        
    end % end of if plotCoG
    
   

end % end of subject Number



beep
