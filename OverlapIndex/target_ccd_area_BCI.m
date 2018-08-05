% This function read the CCD information in the target local area and then find size of the active area.
% function [ActiveRatio,AOI,LNo,AOI_CI,LNo_CI]=target_ccd_area_BL1(subject)
function [OverArea_s, ActiveRatio_R, ActiveRatio_H, maxStrenght_H, maxStrenght_R, CoG_H, CoG_R]=target_ccd_area_BCI(subject)
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
taskList={'HCL-NOT';'HOP-NOT'};
% taskList={'HCL-T';'HOP-T'};
phaseList={'whole'};

plotCort=1;

show_frame=0;
yjcolor=['b','r','m'];
currents={};
V1='BL ';
V2='EXE';
debug_yj=1;
% group={};
timePoints=20;
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
        
        
        for methodInd=1:length(method)
            ActCount=0;
            cdr_file=[cur_task,'.cdr'];
            cdr_BL_file=[cur_task,'_BL','.cdr'];
	%             cdr_file=[cur_task,'.cdr'];
            mean_ccd_file=[Base_dir,subjectName,'\new\',cur_task,'.dat'];
            
            disp('********************');
            disp(cdr_file)
			cdr_file_name=[Base_dir,subjectName,'\new\',cdr_file];
			[cortexL,targetL,targetC]=find_target_PC4(cdr_file_name,mask_file_name,timePoints,plotCort);
            
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
            targetC_H=targetC(:,end-7:end);
            targetC_R=targetC(:,1:13);
            
            [locationNo,TimeNo]=size(targetC);
            
            targetC_H_mean=mean(targetC_H,2);
            maxStrenght_H(i)=max(targetC_H_mean);
            targetC_H_norm=targetC_H_mean/maxStrenght_H(i);
            
            targetC_R_mean=mean(targetC_R,2);
            maxStrenght_R(i)=max(targetC_R_mean);
            targetC_R_norm=targetC_R_mean/maxStrenght_R(i);
            
            save_data=[targetC_R_mean,targetC_H_mean];
            save (mean_ccd_file,'save_data','-ascii');
            
            CoG_H(i,:)=sum(repmat(targetC_H_mean,1,3).*targetL,1)/sum(targetC_H_mean);
            CoG_R(i,:)=sum(repmat(targetC_R_mean,1,3).*targetL,1)/sum(targetC_H_mean);
           
            
            if debug_yj
                figure
                hold on
                h2=plot3(cortexL(:,1),cortexL(:,2),cortexL(:,3),'.','MarkerSize',3);
                CI = griddata(targetL(:,1),targetL(:,2),targetC_H_mean/7,XI,YI);
                surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
                h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
                h1=plot3(CoG_H(i,1),CoG_H(i,2),CoG_H(i,3),'m*','MarkerSize',6);
                view(90,90)
                title (['After the onset of wrist muscles for task: ',cur_task]);
                xlabel('x');ylabel('y');
                colorbar
            end
            
            if debug_yj
                figure
                hold on
                h2=plot3(cortexL(:,1),cortexL(:,2),cortexL(:,3),'.','MarkerSize',3);
                CI = griddata(targetL(:,1),targetL(:,2),targetC_R_mean/14,XI,YI);
                surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
                h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
                h1=plot3(CoG_R(i,1),CoG_R(i,2),CoG_R(i,3),'m*');
                view(90,90)
                title (['Before the onset of wrist muscles for task: ',cur_task]);
                xlabel('x');ylabel('y');
                colorbar
            end
            
            [muhat, muci] = expfit(targetC_H_norm);
            ind=find(targetC_H_norm>muci(1));
            ActiveRatio_H(k,i)=length(ind)/locationNo; %ActiveRatio: the percentage of actived area

            [muhat, muci] = expfit(targetC_R_norm);
            ind=find(targetC_R_norm>muci(1));
            ActiveRatio_R(k,i)=length(ind)/locationNo; %ActiveRatio: the percentage of actived area
           
        end %end of methodInd
        OverArea_norm=targetC_H_norm.*targetC_R_norm; %Overlapped area after normalization of the strength
        [muhat, muci] = expfit(OverArea_norm);
        ind1=find(OverArea_norm>muci(1));
        OverArea_s(k,i)=length(ind1)/locationNo;  % OverArea_s: the overlapped index: overlaped area significantly larger than 0 / total area
    end % end of task
  
    if show_frame
        figure 
        CI = griddata(targetL(:,1),targetL(:,2),OverArea_norm,XI,YI);
        hold on
        surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
        h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
        colorbar
        set(h,'LineWidth',7)
        view(90,90)
    end

    clear OverArea_norm XI YI ZI CI
end % end of subject Number

beep
