% This function read the CCD information in the target local area.
% Input: Subject: This is a colunmm cell. Each element is a string showing
% subject's name
function [OI_time,OI_task]=target_ccd_oi_time_new(subject)
% function target_ccd (subject,cdr_file,edge_file,timePoints,plotCort)
% function plot_target_ccd (subjectName,cdr_file,edge_file,timePoints,plotCort)
% % method={'1','2.5','3','4','5','6','7'};
% % pre_method=method(1);
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
% currents={};
OI_time_allSub=[];
OI_task_allSub=[];
for k=1:subNo
    cur_CI={};
    subjectName=subject{k}
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
    

%     record=[];
    for i=1:length(taskList)
        cur_task=taskList{i};
        
        cdr_file=[cur_task,'.cdr'];
        cdr_BL_file=[cur_task,'_BL','.cdr'];
        disp('********************');
        disp(cdr_file)
		cdr_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_file];
		cdr_BL_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_BL_file];
        
        [Matrix,count,NR]=read_Curry_file3(cdr_file_name,'LOCATION',0,0);
        TimePTS=NR(2);
		[cortexL,targetL,targetC]=find_target_PC3(cdr_file_name,edge,TimePTS,plotCort);
        [Matrix,count,NR]=read_Curry_file3(cdr_BL_file_name,'LOCATION',0,0);
        BL_TimePTS=NR(2);
        [cortexL,targetL,targetC_BL]=find_target_PC3(cdr_BL_file_name,edge,BL_TimePTS,plotCort);
        targetC_BL=mean(targetC_BL,2);
        
        minX=min(targetL(:,1));
        maxX=max(targetL(:,1));
        minY=min(targetL(:,2));
        maxY=max(targetL(:,2));
        minZ=min(targetL(:,3));
        maxZ=max(targetL(:,3));
        x=[minX-5:0.8:maxX+5];
        y=[minY-5:0.8:maxY+5];
        [XI,YI] = meshgrid(x,y); 
        

        CI_BL = griddata(targetL(:,1),targetL(:,2),targetC_BL,XI,YI);
		TA=isnan(CI_BL);
		ind=find(TA==1);
		CI_BL(ind)=0;
        
        winOverlap=5;
        winSize=10;
        LoopNo=ceil((TimePTS-winOverlap)/(winSize-winOverlap));
        for phTime=1:LoopNo
            ZI = griddata(targetL(:,1),targetL(:,2),targetL(:,3),XI,YI);
            if show_frame figure (i); end
            record=[];
            position=[];
            k=1;
            sum_CI=zeros(size(XI));
            for timenow=(phTime-1)*(winSize-winOverlap)+1:(phTime-1)*(winSize-winOverlap)+winSize
                CI = griddata(targetL(:,1),targetL(:,2),targetC(:,timenow),XI,YI);
				TA=isnan(CI);
				ind=find(TA==1);
				CI(ind)=0;
%                 CI=CI-CI_BL;

                sum_CI=sum_CI+CI;
                
                strength_now=targetC(:,timenow);
                strength_norm=strength_now/max(strength_now);
                ind=find(strength_norm<0.25);
                strength_now(ind)=0;
                
                center_now=sum(repmat(strength_now,1,3).*targetL,1);
                center_now=center_now/sum(strength_now);
                center_str_now=mean(strength_now);
                
                [dis2censul,index]=min(sqrt( sum( (censul-repmat(center_now,size(censul,1),1)).^2,2 )  ));
                if dis2censul<200
                    midL=sqrt( sum( (censul(1,:)-center_now).^2,2 )  );
                    if center_now(2)<censul(index,2)
                        position(k)=-1;
                    else
                        position(k)=1;
                    end
                    k=k+1;
                    record=[record;(timenow-1)*3.9-100,midL,center_str_now];
                end
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
                end %end of show_frame
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
            if phTime==1
                pre_current=sum_CI./winSize;
            else
                cur_current=sum_CI./winSize;
                OI_time(i,phTime-1)=sum(sum(pre_current.*cur_current,1),2)/sum(sum(pre_current,1),2); %OI stores the overlap between current cortial currents and currents in previous window
                cur_current=pre_current;
            end
            cur_CI{i,phTime}=sum_CI./winSize; %cur_CI stores the currents in current window
        end %end of phTime
    end % end of task
    for phTime=1:LoopNo
        C1=cur_CI{1,phTime};
        C2=cur_CI{2,phTime};
        OI_task(1,phTime)=sum(sum(C1.*C2,1),2)/sum(sum((C2),1),2);
        OI_task(2,phTime)=sum(sum(C1.*C2,1),2)/sum(sum((C1),1),2);
    end
	beep
    OI_time_allSub=[OI_time_allSub,OI_time'];
    OI_task_allSub=[OI_task_allSub,OI_task'];
    outputfile1=['C:\Documents and Settings\Jun Yao\My Documents\Experiment\sh&Elb\timeRes\OI_time.txt'];
    save (outputfile1,'OI_time_allSub','-ascii')
    outputfile2=['C:\Documents and Settings\Jun Yao\My Documents\Experiment\sh&Elb\timeRes\OI_task.txt'];
    save (outputfile2,'OI_task_allSub','-ascii')

end % end of subject Number
% rmpath('C:\Documents and Settings\Jun Yao\Matlab\Eeg\inv')
% figure (3)
% xlabel('time (ms)')
% ylabel('distance (mm)')
% zlabel('strength (uA)')
% legend('SABD','EF')
% title(['threshold=15,0.5',method])
% % legend('ind','sld')
% view(0,90)
% grid on
% grid minor
% % minrow=min([size(r1,1),size(r2,1),size(r3,1)])
% % p=anova1([r3(1:minrow,2),r2(1:minrow,2),r1(1:minrow,2)])
% % p1=anova1([r3(1:minrow,2),r2(1:minrow,2)])
% % p2=anova1([r3(1:minrow,2),r1(1:minrow,2)])
% % p3=anova1([r2(1:minrow,2),r1(1:minrow,2)])
% minrow=min([size(r1,1),size(r2,1)])
% p3=anova1([r2(1:minrow,2),r1(1:minrow,2)])
% ind_shld=[ mean(r1(1:minrow,2)),std(r1(1:minrow,2))/sqrt(timePoints),mean(r2(1:minrow,2)),std(r2(1:minrow,2))/sqrt(timePoints) ]
% % ind_shld=[ mean(r1(1:minrow,2)),std(r1(1:minrow,2))/sqrt(5),mean(r2(1:minrow,2)),std(r2(1:minrow,2))/sqrt(5) ]
% 
% % ind_wst_shld=[ mean(r1(1:minrow,2)),std(r1(1:minrow,2)),mean(r2(1:minrow,2)),std(r2(1:minrow,2)),mean(r3(1:minrow,2)),std(r3(1:minrow,2)), p, p1, p2,p3]
% % legend('ind','wst','sld')
% figure(5)
% xlabel('time (ms)')
% ratio
beep
% legend('ind','wst','sld')
% figure(2)
% title(['threshould=15,0.5 ',method])
% xlabel('sld   wst   ind')
