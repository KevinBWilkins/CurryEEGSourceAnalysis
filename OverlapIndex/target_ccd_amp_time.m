% This function read the CCD information in the target local area.
function [max_amp,avg_amp]=target_ccd_amp_time(subject)
% function target_ccd (subject,cdr_file,edge_file,timePoints,plotCort)
% function plot_target_ccd (subjectName,cdr_file,edge_file,timePoints,plotCort)
method={'1','2.5','3','4','5','6','7'};
pre_method=method(1);
addpath('C:\Documents and Settings\Jun Yao\Matlab\Eeg\inv')

% subject={'brad';'Jules';'TK';'Ty';'LK'};
% subject={'BH';'JD';'Ty';'LK'};
% subject={'JD'};

Base_dir='C:\Documents and Settings\Jun Yao\Subject\';
subNo=length(subject);
% taskList={'ext';'flex';'sabd';'sadd'};
% taskList={'ext';'flex';'abd';'add'};
taskList={'abd';'flex'};
% taskList={'ind';'sld'};
% method='LRT15';
% phaseList=[-700, -250, 600;  -700, -250, 600; -700, -250, 600; -700, -150, 700; -700, -200, 600];
%phaseList={'pre'; 'mov'; 'hd'};  
phaseList={'whole'};

% phaseList=[ -700, -200, 600];
plotCort=0;

show_frame=0;
yjcolor=['b','r','m'];
currents={};
for k=1:subNo
    subjectName=subject{k};
	edge1_file_name=[Base_dir,subjectName,'\Data\MRI\','SM1_L.sp'];
    edge2_file_name=[Base_dir,subjectName,'\Data\MRI\','PREMOTOR_L.sp'];
	edge3_file_name=[Base_dir,subjectName,'\Data\MRI\','SM1_R.sp'];
    edge4_file_name=[Base_dir,subjectName,'\Data\MRI\','PREMOTOR_R.sp'];
	
    [edge1,Ecount,ENR]=read_Curry_file3(edge1_file_name,'LOCATION',0,0);
	[edge2,Ecount,ENR]=read_Curry_file3(edge2_file_name,'LOCATION',0,0);
	[edge3,Ecount,ENR]=read_Curry_file3(edge3_file_name,'LOCATION',0,0);
	[edge4,Ecount,ENR]=read_Curry_file3(edge4_file_name,'LOCATION',0,0);
	edge=[edge1;edge2;edge3;edge4];
    
    censul_file_name=[Base_dir,subjectName,'\Data\MRI\','CS_L.sp'];
    [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name,'LOCATION',0,0);
    censul_file_name=[Base_dir,subjectName,'\Data\MRI\','CS_R.sp'];
    [censul2,Ecount,ENR]=read_Curry_file3(censul_file_name,'LOCATION',0,0);
    censul=[censul1;censul2];
    record=[];
    for i=1:length(taskList)
        cur_task=taskList{i};
        for phTime=1:length(method)
            timePoints=13;
            cdr_file=[cur_task,'_',method{phTime},'.cdr'];
            disp('********************');
            disp(cdr_file)
			cdr_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_file];
			
			[cortexL,targetL,targetC]=find_target_PC3(cdr_file_name,edge,timePoints,plotCort);

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
            if show_frame figure (i); end

%             Tes = delaunay3(targetL(:,1),targetL(:,2),targetL(:,3));
            record=[];
            position=[];
            k=1;
            sum_CI=zeros(size(XI));
            max_C=0;
            for timenow=1:timePoints
                CI = griddata(targetL(:,1),targetL(:,2),targetC(:,timenow),XI,YI);
				TA=isnan(CI);
				ind=find(TA==1);
				CI(ind)=0;
                [CI_xDim,CI_yDim]=size(CI);

                sum_CI=sum_CI+CI;
                max_C=max_C+max(targetC(:,timenow));
%                 strength_now=targetC(:,timenow);
%                 strength_norm=strength_now/max(strength_now);
%                 ind=find(strength_norm<0.25);
%                 strength_now(ind)=0;
%                 
%                 center_now=sum(repmat(strength_now,1,3).*targetL,1);
%                 center_now=center_now/sum(strength_now);
%                 center_str_now=mean(strength_now);
%                 
%                 [dis2censul,index]=min(sqrt( sum( (censul-repmat(center_now,size(censul,1),1)).^2,2 )  ));
%                 if dis2censul<200
%                     midL=sqrt( sum( (censul(1,:)-center_now).^2,2 )  );
%                     if center_now(2)<censul(index,2)
%                         position(k)=-1;
%                     else
%                         position(k)=1;
%                     end
%                     k=k+1;
%                     record=[record;(timenow-1)*3.9-100,midL,center_str_now];
%                 end
%                 if show_frame
%                    C=targetC(:,timenow)/max(targetC(:,timenow));
% %                    subplot(4,5,timenow)
% 	%                 plot3(cortexL(:,1),cortexL(:,2),cortexL(:,3),'b.');
%                     hold on
% 	%                 patch('Vertices',targetL,'Faces',Tes,'FaceVertexCData',targetC(:,timenow)/max(targetC(:,timenow)),'FaceColor','interp','EdgeColor','none','FaceAlpha',0.2);
%                     CI = griddata(targetL(:,1),targetL(:,2),C,XI,YI);
%                     surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6)
% %                     plot3(censul(:,1),censul(:,2),censul(:,3),'LineWidth',7,'color','g');
%                     h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
%                     set(h,'LineWidth',7)
%                     plot3(center_now(:,1),center_now(:,2),center_now(:,3),'m*')
%                     view(90,90)
%                     
% %                     view(-90,90)
%                     startTime=-100*str2num(method);
%                     title([cur_task,'    Time=',num2str(startTime),'(ms)'],'FontSize',8); 
% 	%                 pause
% 	%                 clf
%                 end %end of show_frame
            end% end of timenow
%             switch i
%                 case 1
%                     plot3(record(:,1),record(:,2),record(:,3),'b.','MarkerSize',20)
%                     r1=record;
%                 case 2
%                     plot3(record(:,1),record(:,2),record(:,3),'r.','MarkerSize',20)
%                     r2=record;
%                 case 3
%                     plot3(record(:,1),record(:,2),record(:,3),'m.','MarkerSize',20)
%                     r3=record;
%             end
%             if phTime==1
%                 pre_current=sum_CI./timePoints;
%             else
%                 cur_current=sum_CI./timePoints;
%                 OI(i,phTime-1)=sum(sum(pre_current.*cur_current,1),2)/sum(sum(pre_current,1),2);
%                 cur_current=pre_current;
%             end
            avg_amp(i,phTime)=sum(sum(sum_CI,1),2)/(CI_xDim*CI_yDim-length(ind))/timePoints;
            max_amp(i,phTime)=max_C/timePoints;

        end %end of phTime
    end % end of task
end % end of subject Number
beep
