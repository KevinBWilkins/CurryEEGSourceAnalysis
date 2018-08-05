% This function read the CCD information in the target local area.
function AOI=target_ccd_cog_aoi1(subject)
% function target_ccd (subject,cdr_file,edge_file,timePoints,plotCort)
% function plot_target_ccd (subjectName,cdr_file,edge_file,timePoints,plotCort)

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
    
    censul_file_name=[Base_dir,subjectName,'\Data\MRI\','CS_L.sp'];
    [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name,'LOCATION',0,0);
    censul_file_name=[Base_dir,subjectName,'\Data\MRI\','CS_R.sp'];
    [censul2,Ecount,ENR]=read_Curry_file3(censul_file_name,'LOCATION',0,0);
    censul=[censul1;censul2];
    record=[];
    for i=1:length(taskList)
        cur_task=taskList{i};
        for phTime=1:length(phaseList)
            timePoints=155;
            cdr_file=[cur_task,'.cdr'];
            disp('********************');
            disp(cdr_file)
			cdr_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_file];
			
% 			[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
			[cortexL,targetL,targetC]=find_target_PC1(cdr_file_name,edge1_file_name,edge2_file_name,edge3_file_name,edge4_file_name,timePoints,plotCort);

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

            record=[];
            position=[];
            sum_CI=zeros(size(XI));
            winPoints=13;
            constant_time=155-1*winPoints;
            constant_time=0;

            for timenow=1:winPoints
                C=targetC(:,timenow+constant_time)/max(targetC(:,timenow+constant_time));
                CI = griddata(targetL(:,1),targetL(:,2),C,XI,YI);
				TA=isnan(CI);
				ind=find(TA==1);
				CI(ind)=0;

                sum_CI=sum_CI+CI;
                
                strength_now=targetC(:,timenow+constant_time);
                strength_norm=strength_now/max(strength_now);
                ind=find(strength_norm<0.25);
                strength_now(ind)=0;
                
                center_now=sum(repmat(strength_now,1,3).*targetL,1);
                center_now=center_now/sum(strength_now);
                center_str_now=mean(strength_now);
                [dis2censul,index]=min(sqrt( sum( (censul-repmat(center_now,size(censul,1),1)).^2,2 )  ));
                if show_frame
                   C=targetC(:,timenow)/max(targetC(:,timenow));
%                    subplot(4,5,timenow)
	%                 plot3(cortexL(:,1),cortexL(:,2),cortexL(:,3),'b.');
                    hold on
	%                 patch('Vertices',targetL,'Faces',Tes,'FaceVertexCData',targetC(:,timenow)/max(targetC(:,timenow)),'FaceColor','interp','EdgeColor','none','FaceAlpha',0.2);
                    CI = griddata(targetL(:,1),targetL(:,2),C,XI,YI);
                    surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6)
%                     plot3(censul(:,1),censul(:,2),censul(:,3),'LineWidth',7,'color','g');
                    h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
                    set(h,'LineWidth',7)
                    plot3(center_now(:,1),center_now(:,2),center_now(:,3),'m*')
                    view(90,90)
                    
%                     view(-90,90)
                    startTime=-100*str2num(method);
                    title([cur_task,'    Time=',num2str(startTime),'(ms)'],'FontSize',8); 
	%                 pause
	%                 clf
                end
            end
        end %end of time
        current{i}=sum_CI./winPoints;
    end % end of task
    C1=current{1};
    C2=current{2};
    AOI(k)=sum(sum(C1.*C2,1),2)/size(C2,1)/size(C2,2);
end % end of subject Number
