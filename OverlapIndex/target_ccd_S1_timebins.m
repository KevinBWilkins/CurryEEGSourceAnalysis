% This function read the CCD information in the target local area during a time window and then find size of the active area.
% function [ActiveRatio,AOI,LNo,AOI_CI,LNo_CI]=target_ccd_area_BL1(subject)
function [act_ratio_loc_time, LI]=target_ccd_S1_timebins(subject,group,arm)
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
% taskList={'shoulder-l-11';'midfinger-l-5'};%}; %};%;'Lip-2'};%'Midfinger-l-7'; %Tim pre-TR
% taskList={'shoulder-l-9';'mid-fig-3';'shoulder-r-12';'mid-r-4'};%'Lip-3'}; %tim post1 
% taskList = {'Lsh-9';'needle1';'Rsh-9';'RmidFinger'};
% taskList={'sup_Lshoulder';'sup_Lelbow';'sup_Rshoulder';'sup_Rfing';};
% taskList={'Ne_Lelbow';'Sup_Llfinger';'Ne_Relbow';'Sup_Rlfinger'};
% taskList={'Ne_Lelbow2';'Sup_Lelbow2';'Sup_Llfinger2';'Ne_Relbow2';'Sup_Re
% lbow2';'Sup_Rlfinger2'};
taskList={'Sup_Rlfinger2-2'};
% taskList={'Ne_Lelbow';'Cu_Llfinger';'Cu_Llfinger2';'Ne_Relbow';'Cu_Rlfinger';'Cu_Rlfinger2'};
% taskList={'Sup_Llfinger2'};
% taskList={'Ne_Lelbow';'Ne_Relbow';'Sup_Lelbow';'Sup_Lfing';'Sup_Relbow'};%'Lsh-9';'needle1';'rbic-4';'RmidFinger';'Rsh-9'};
% taskList={'back_hand';'Ind+thumb_new';'L-shoulder';'R-back_hand';'R-index';'R-shoulder_new';'R-index-mirrow_new'};
% taskList={'Ind+thumb_new'};

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


winOverlap=0;
winSize=1;
timePoints=65;
POINT_T_DELTA        =  0.49;
start_T = [5.97,6.46,6.95,6.95,6.95,5.97];
start_T = [9.89];
piont_shift = round((10-start_T)/POINT_T_DELTA);


LoopNo=fix((timePoints-winOverlap)/(winSize-winOverlap));
needed_len=(LoopNo-1)*(winSize-winOverlap)+winSize;
startPoints=timePoints-needed_len-1;
taskColorMap = ['jet';'jet'];
for k=1:subNo
%     subjectArm=arm{k}; 
%     if subjectArm == 'l'
%         location_order = [3 1 7 5 6 8 2 4];
%     elseif subjectArm == 'r'
        location_order = [6 8 2 4 3 1 7 5];   
%     end
    
    subjectName=subject{k};   
    %--------------------------------------------
    %Read in central sulcus file, may be different name
    censul_file_name1=[Base_dir,subjectName,'\Points\','CS_cc.pom'];   %central sulcus
    mask_file_name=[Base_dir,subjectName,'\Points\','Mask'];   %index of ROI
%     mask_file_name=[Base_dir,subjectName,'\Points\','Mask-pre'];   %index of ROI

    
    

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
        cdr_file=[cur_task,'.cdr'];
        cdr_BL_file=[cur_task,'_BL','.cdr'];
%             cdr_file=[cur_task,'.cdr'];
        mean_ccd_file=[Base_dir,subjectName,'\new\',cur_task,'.dat'];
        
        disp('********************');
        disp(cdr_file)
        cdr_file_name=[Base_dir,subjectName,'\new\',cdr_file];
%         cdr_file_name=[Base_dir,subjectName,'\AKpre\results_3mm\',cdr_file];
% 			[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
%         [cortexL,targetL_all,targetC_wh_time_ROI(:,i)]=find_target_PC4(cdr_file_name,mask_file_name,timePoints,plotCort,0:1);
        [cortexL,targetL_all,targetC_wh_time_ROI(:,i)]=find_target_PC4(cdr_file_name,mask_file_name,timePoints+piont_shift(i),plotCort,0:1); %this is only for Tim post3

        targetC_all = cell2mat(squeeze(targetC_wh_time_ROI(:,i)));
%         targetC_all = targetC_all (:, [21:end]);
        max_targetC(i) = max(max(targetC_all));
        all_targetC(i) = sum (sum(targetC_all));
        total_location_no = size(targetC_all,1);
%         load targetC_wh_time
        if (i==1 )
                all_targetL = cell2mat(targetL_all);
%                 all_targetL = targetL_all_mat;
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
%-----separate location in S1 and change targetL_all to cell_array  
    
    
%----------------------------------------------------------------------    
%     for phTime=1:LoopNo+1
        
%         if phTime <= LoopNo
%             cur_window = (phTime-1)*(winSize-winOverlap)+1+startPoints:(phTime-1)*(winSize-winOverlap)+winSize+startPoints;
%         else
%             cur_window = (phTime-1)*(winSize-winOverlap)+1+startPoints:timePoints;
%            
%         end
%         cur_window = phTime;
        for i=1:length(taskList)
            cur_task=taskList{i}
            [targetL_all_div,cur_targetC_wh_time_ROI_div] = divide_loc_current(targetL_all,squeeze(targetC_wh_time_ROI(:,i)), 1:2, piont_shift(i));
            [m,n]=size(targetL_all_div);
            loc_div =0;
            for i_location=1:m % location
                 for i_div = 1 : n
                    if i_location == 1
                        div_ind = n-i_div+1;
                    else
                        div_ind = i_div;
                    end

                    targetC=cur_targetC_wh_time_ROI_div{i_location,div_ind};
%                     targetC=targetC_cell(:,cur_window);
                
%                     targetL=targetL_all_div{i_location,div_ind};

                    [locationNo]=size(targetC,1);
                    loc_div = loc_div+1;
%---------------------------------calculate the acitvation ratio-----------------
                    if locationNo
%                         act_ratio_loc_time (k,i,loc_div,:)= sum(targetC,1)*total_location_no/all_targetC(i)/locationNo;
                        act_ratio_loc_time (k,i,loc_div,:)= sum(targetC,1)/all_targetC(i)/locationNo;
                        
                    else
                        act_ratio_loc_time (k,i,loc_div,:)= zeros(1,timePoints);
                    end
%______________________________________________________________________________________                
%                 if debug_yj
%                     if cur_task(1)=='a'
%                         figure (2) %*phTime)
% 
%                     else
%                         figure (2*phTime+1)
%                     end
%                     clf
%                     hold on
%                     CI = griddata(targetL(:,1),targetL(:,2),targetC_mean,XI,YI);
%                     surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
%                     h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
%                     h1=plot3(CoG(k,phTime,i,1),CoG(k,phTime,i,2),CoG(k,phTime,i,3),'m*','MarkerSize',6);
%                     view(90,90)
%                     colorbar
%                     h=plot3(activeL(:,1),activeL(:,2),activeL(:,3),'mo');
%                     title (cur_task);
% 
%                 end
%                 if phTime~=1 fig_h=plot_data_on_cortex(subjectName, -21, Base_dir, targetL, censul, targetC_norm, cortexL, phTime-1, taskColorMap(i,:), i, cur_window);  end
                end %end of division index
            end %end of the location
            %------------------the following part only fit to Tim's results
            %with the 2 left sites first then followed by the 2 right sites
            if i<0 %left side
                tmp_i = squeeze(act_ratio_loc_time(k,i,[1:9],:));
                tmp_c = squeeze(act_ratio_loc_time(k,i,[10:18],:));
            else i>0 %right side
                tmp_c = squeeze(act_ratio_loc_time(k,i,[1:9],:));
                tmp_i = squeeze(act_ratio_loc_time(k,i,[10:18],:));
            end
            %--------------------------------------------------------------
            
            LI(k,i,:)= (sum(tmp_c)-sum(tmp_i)) ./ (sum(tmp_c)+sum(tmp_i));
        end % end of task

    clear XI YI ZI CI
    X=0:POINT_T_DELTA:(timePoints-1)*POINT_T_DELTA;
    X=X+10;
    X=fliplr(X);

    plot_act_time_loc_S1 (act_ratio_loc_time, LI);


    output_file_name=[Base_dir,subjectName];%, taskList{i}];        

    save (output_file_name,'act_ratio_loc_time','LI')

    
   

end % end of subject Number


