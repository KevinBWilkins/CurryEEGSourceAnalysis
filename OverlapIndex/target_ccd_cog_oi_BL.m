% This function read the CCD information in the target local area.
function target_ccd_cog_oi_BL(method)
% function target_ccd (subject,cdr_file,edge_file,timePoints,plotCort)
% function plot_target_ccd (subjectName,cdr_file,edge_file,timePoints,plotCort)

addpath('C:\Documents and Settings\Jun Yao\Matlab\Eeg\inv')

% subject={'brad';'Jules';'TK';'Ty';'LK'};
% subject={'BH';'JD';'Ty';'LK'};
% subject={'MG';'PG';'LN';'VP';'JW';'FN';'FRG'};
% subject={'Dan';'MLM';'LZ';'JD';'TK';'BH';'BS'};
subject={'FN2'};
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
    record=[];
    for i=1:length(taskList)
        cur_task=taskList{i};
        for phTime=1:length(phaseList)
%             if (k<3 & phTime==2)
%                 timePoints=65;
%             else
%                 timePoints=77;
%             end
            timePoints=155;
%             cdr_file=[subjectName,'_',cur_task,'_',num2str(phaseList(k,phTime)),'.cdr'];
%             cdr_file=[subjectName,cur_task,'_',phaseList{phTime},'.cdr'];
%             cdr_file=[cur_task,'_',method,'.cdr'];
%             cdr_file=[cur_task,'_',method,'.cdr'];
            cdr_file=[cur_task,'.cdr'];
            cdr_BL_file=[cur_task,'_BL','.cdr'];


            disp('********************');
            disp(cdr_file)
			cdr_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_file];
			cdr_BL_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_BL_file];

% 			[cortexL,targetL,targetC]=find_target_PC(cdr_file_name,edge1_file_name,edge2_file_name,timePoints,plotCort);
% 			[cortexL,targetL,targetC]=find_target_PC1(cdr_file_name,edge1_file_name,edge2_file_name,edge3_file_name,edge4_file_name,timePoints,plotCort);
            [cortexL,targetL,targetC]=find_target_PC3(cdr_file_name,edge,timePoints,plotCort);
            [cortexL,targetL,targetC_BL]=find_target_PC3(cdr_BL_file_name,edge,13,plotCort);
            targetC_BL=mean(targetC_BL,2);

            targetC=targetC-repmat(targetC_BL,1,timePoints);
            
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
            kk=1;
            sum_CI=zeros(size(XI));
            winPoints=13;
            constant_time=155-2*winPoints;
            
            for timenow=1:winPoints
%             for timenow=9:9
                
                CI = griddata(targetL(:,1),targetL(:,2),targetC(:,timenow+constant_time),XI,YI);
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
%                 [center_str_now,index]=max(strength_now);
%                 if length(index)>1
%                     mindis=100;
%                     for k=1:length(index)
%                         dis=min(sqrt( sum( (censul-repmat(target(index,:),size(censul,1),1)).^2,2 )  ));
%                         if dis<mindis
%                             ind_index=k;
%                             mindis=dis;
%                         end
%                     end
%                     index1=index(ind_index);
%                 else
%                     index1=index;
%                 end
%                 center_now=targetL(index1,:);
                
                [dis2censul,index]=min(sqrt( sum( (censul-repmat(center_now,size(censul,1),1)).^2,2 )  ));
                if dis2censul<200
                    midL=sqrt( sum( (censul(1,:)-center_now).^2,2 )  );
                    if center_now(2)<censul(index,2)
                        position(kk)=-1;
                    else
                        position(kk)=1;
                    end
                    kk=kk+1;
                    record=[record;(timenow-1)*3.9-100,midL,center_str_now];
                end
                if show_frame
                   C=targetC(:,timenow+constant_time)/max(targetC(:,timenow+constant_time));
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
%                     view(90,90)
                    
                    view(-90,90)
                    startTime=-100*str2num(method);
                    title([cur_task,'    Time=',num2str(startTime),'(ms)'],'FontSize',8); 
	%                 pause
	%                 clf
                end
            end
%             t=10:0.488:(timePoints-1)*0.488+10;
%             midL=0:round(sqrt( sum( (censul(:,end)-censul(:,1)).^2,2 ))+1);
%             [TI,midLI]=meshgrid(t,midL);
%             result = griddata(record(:,1),record(:,2),record(:,3),TI,midLI);
%             mesh(t,midL,result);
%             title(cur_task)
%             xlabel('time (ms)')
%             ylabel('distance (mm)')
%             zlabel('strength (uA)')
            targetC_file=[Base_dir,subjectName,'\results_inv\noVMOV\',subjectName,'_',cur_task,'_',method,'.mat']
            save (targetC_file,'record');
            figure(3)
            hold on
            
            switch i
                case 1
                    plot3(record(:,1),record(:,2),record(:,3),'b.','MarkerSize',20)
                    r1=record;
                case 2
                    plot3(record(:,1),record(:,2),record(:,3),'r.','MarkerSize',20)
                    r2=record;
                case 3
                    plot3(record(:,1),record(:,2),record(:,3),'m.','MarkerSize',20)
                    r3=record;
            end
            figure(5)
            hold on
            subplot(length(taskList),1,i)
            plot(record(:,1),position,yjcolor(i));
            axis([10,40,-1.5,2])
            legend(taskList{i})
            positive=find(position>0);
			ratio(i)=length(positive)/length(position);
            if i==1
                title(['Position (post/pre-central gyrus threshould=0.25 ',method,')'])
            end
            ylabel('position')
            text(12,1.5,['ratio=',num2str(ratio(i))],'FontSize',8);
        end %end of time
        current{i}=sum_CI./winPoints;
%         current{i}=sum_CI./5;
    end % end of task
    C1=current{1};
    C2=current{2};
    OI1_2(k)=sum(sum(C1.*C2,1),2)/sum(sum(C2,1),2)
    OI2_1(k)=sum(sum(C1.*C2,1),2)/sum(sum(C1,1),2)
    
	beep
end % end of subject Number
% rmpath('C:\Documents and Settings\Jun Yao\Matlab\Eeg\inv')
figure (3)
xlabel('time (ms)')
ylabel('distance (mm)')
zlabel('strength (uA)')
legend('SABD','EF')
title(['threshold=15,0.5',method])
% legend('ind','sld')
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
% beep
% legend('ind','wst','sld')
% figure(2)
% title(['threshould=15,0.5 ',method])
% xlabel('sld   wst   ind')
