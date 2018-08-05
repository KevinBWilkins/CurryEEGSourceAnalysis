% This function read the CCD information in the target local area and then find the CoG of that part.
% function [ActiveRatio,AOI,LNo,AOI_CI,LNo_CI]=target_ccd_area_BL1(subject)
function [CoG_s]=target_ccd_onetime(subject,group,task,startTime)

method={'1'};
% pre_method=method(1);
addpath('C:\Documents and Settings\Jun Yao\Matlab\Eeg\inv')


Base_dir=['F:\data\inverse_results\',group,'\'];
subNo=length(subject);
% phaseList={'AKpre_CC';'post'};
% phaseList={'AKpre_CC'};
phaseList={'post'};


plotCort=1;

show_frame=0;
yjcolor=['b','r','m'];
currents={};
V1='BL ';
V2='EXE';
debug_yj=1;
% group={};
save_data=[];
alpha_static = 0.01;

for k=1:subNo
    subjectName=subject{k};   
    %--------------------------------------------
    %Read in central sulcus file, may be different name
    censul_file_name1=[Base_dir,subjectName,'\Points\','CS_cc.pom'];   %central sulcus
%     mask_file_name=[Base_dir,subjectName,'\Points\','Mask-l-smc.mat'];   %index of ROI
    mask_file_name=[Base_dir,subjectName,'\Points\','Mask-sensorimotor.mat'];   %index of ROI
    
    %if only one central sulcus file
    [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name1,'LOCATION',0,0);
    censul2 = censul1;
    censul=[censul1];
    %----------------------------------------------------------------------

    record=[];
    winPoints=15;
    for i=1:length(phaseList)
        
        cur_phase=phaseList{i};
        timePoints=[31];
        
        ActCount=0;
        cdr_file=[task,'.cdr'];
%         cdr_BL_file=[cur_task,'_BL','.cdr'];
%             cdr_file=[cur_task,'.cdr'];
        mean_ccd_file=[Base_dir,subjectName,'\',cur_phase,'\',task,'_',num2str(startTime),'_onetime.dat'];
        
        disp('********************');
        disp(cdr_file)
        cdr_file_name=[Base_dir,subjectName,'\',cur_phase,'\',cdr_file];
%         cdr_BL_file_name=[Base_dir,subjectName,'\new\',cdr_BL_file];
% 			[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
        [cortexL,targetL,targetC]=find_target_PC5(cdr_file_name,mask_file_name,timePoints (i),plotCort,0);

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
        targetC=targetC(:,startTime(k,i));
%         targetC=targetC(:,startTime(k,i):startTime(k,i)+winPoints);
%         targetC_BL=[];
%         [cortexL,targetL,targetC_BL]=find_target_PC4(cdr_BL_file_name,mask_file_name,25,plotCort);
%         targetC_BL=mean(targetC_BL,2);
% 			[cortexL,targetC]=find_strength_whole(cdr_file_name,timePoints);
% 			[cortexL,targetC_BL]=find_strength_whole(cdr_BL_file_name,timePoints);

        [locationNo,TimeNo]=size(targetC);

%         targetC=targetC-repmat(targetC_BL,1,winPoints);
        targetC_mean=mean(targetC,2);
        
%----------save the mean for showing purpose---------------
        save_data=[save_data,targetC_mean];
        save (mean_ccd_file,'save_data','-ascii');
%---------------------------------------------------------------
        
        
%         maxStrenght(k,i)=max(targetC_mean);
%         targetC_norm=targetC_mean/maxStrenght(k,i);
% %         targetC_BL_norm=targetC_BL/maxStrenght(k,i);
%         
% 
%         [muhat, muci] = expfit(targetC_norm,alpha_static)
%         thresh=muci(2);
% %         thresh=3*muci(2);
% %         thresh=0.3*max(targetC_norm);
%         
%         ind = find(targetC_norm>thresh);
%         active_ind(i,:) = (i+1)*ones (locationNo, 1);
%         active_ind(i,ind)= 0;
%         
% %         targetC_norm(ind)=0;
% %         targetC_mean(ind)=0;
%         ActiveRatio(k,i)=length(ind)/locationNo; %ActiveRatio: the percentage of actived area
%         
%         activeC = targetC_mean (ind,:);
%         activeL = targetL (ind,:);
        CoG(k,i,:)=sum(repmat(targetC_mean,1,3).*targetL,1)/sum(targetC_mean);
%         CCD_task(:,i)=targetC_mean;
%         CCD_task_norm(:,i)=targetC_norm;
        
        if debug_yj
            if cur_phase(4)=='r'
                figure (2)
               
            else
                figure (4)
            end
            clf
            hold on
            CI = griddata(targetL(:,1),targetL(:,2),targetC_mean,XI,YI);
            surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
            h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
            h1=plot3(CoG(k,i,1),CoG(k,i,2),CoG(k,i,3),'m*','MarkerSize',6);
            view(90,90)
            colorbar
%             h=plot3(activeL(:,1),activeL(:,2),activeL(:,3),'mo');
            title ([cur_phase,'_',task]);
%             subject_pic_Dir=[Base_dir,subjectName,'\pic'];
%             if ~exist(subject_pic_Dir,'dir')
%                 mkdir(subject_pic_Dir);
%             end
%             fig_file_name=[Base_dir,subjectName,'\pic\',cur_task,'.jpg'];   %index of ROI
%             saveas (h, fig_file_name, 'jpg');
            
        end

    end % end of task
%     OverArea_norm=CCD_task_norm(:,1).*CCD_task_norm(:,2); %Overlapped area after normalization of the strength
%     OverArea=CCD_task(:,1).*CCD_task(:,2); %Overlapped area without normalization
%     OI(k,1)=sum(OverArea)/sum(CCD_task(:,2)); %Overlap Index: Overlapped area / the total area of the secondary direction
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
%     OI_inter(k,1)=sum(sum(OverArea_inter,1),2)/sum(sum(CI_task2,1),2); %Overlap index: Overlapped area / the total area of the secondary direction; Differently, first interplate data to a grid
%     OI_inter(k,2)=sum(sum(OverArea_inter,1),2)/sum(sum(CI_task1,1),2);

%     [muhat, muci] = expfit(OverArea_norm, alpha_static);
%     thresh=muci(2);
% %     thresh=3*muci(2);
% %     thresh=0.3*max(OverArea_norm);
%     ind3=find(OverArea_norm>thresh);
%     OverArea_r(k)=length(ind3)/locationNo;  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
%     OAA(k)=2*length(ind3)/(sum(ActiveRatio,2)*locationNo);  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
%     
%     diff=active_ind(1,:)-active_ind(2,:);
%     ind4=find(diff==0);
%     OverArea_r1(k)=length(ind4)/locationNo;  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
%     OAA1(k)=2*length(ind4)/(sum(ActiveRatio,2)*locationNo);  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
%     
%     
%     if debug_yj
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
% 
%     end
%     
%     if show_frame
%         figure 
%         CI = griddata(targetL(:,1),targetL(:,2),OverArea,XI,YI);
%         hold on
%         surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
%         h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
%         colorbar
%         set(h,'LineWidth',7)
%         view(90,90)
%     end
%     CoG_s=['(',num2str(CoG(:,1,1)), ', ',num2str(CoG(:,1,2)),', ',num2str(CoG(:,1,3)),')';'(',num2str(CoG(:,2,1)),', ',num2str(CoG(:,2,2)),', ',num2str(CoG(:,2,3)),')'];
    CoG_s={num2str(CoG(:,1,1)),num2str(CoG(:,1,2)),num2str(CoG(:,1,3));num2str(CoG(:,2,1)),num2str(CoG(:,2,2)),num2str(CoG(:,2,3))};
    clear CCD_task CCD_task_norm OverArea OverArea_norm XI YI ZI CI
    
end % end of subject Number

beep
