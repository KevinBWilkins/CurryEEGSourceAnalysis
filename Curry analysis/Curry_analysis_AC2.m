% This function read the CCD information in the target local area and then find size of the active area.
% function [ActiveRatio,OI,OI_interp,AOI_L]=target_ccd_area_BL1(subject)
% 
% OI is the overlap index taking accout both the area and amplitude
% ActiveRatio is ratio between active location and the total SM area
% OverArea_s is the overlaped active area (area only), no threshould is
%          applied
% AOI_L is the overlaped active location where the activity
%          is larger than 0.25
% AOI_inter is the overlaped area (i.e., obtained based on the interplation
%          on the ective location.) where the activity is larger than 0.25
%
% INSTRUCTIONS ON HOW TO USE:
%
% 1. Starting the function- there are 2 options
%   a. Run the program from the matlab editor
%   b. Run the program from the command line:
%       i. [OI,OI1,ActiveRatio,OverArea_s,AOI_L,AOI_inter,Subjects]=target_ccd_area_BL2()
%           This basically runs the program without arguments. You will
%           need to enter what the subject's name is within the function
%           itself.
%       ii. [OI,OI1,ActiveRatio,OverArea_s,AOI_L,AOI_inter,Subjects]=target_ccd_area_BL2(Subjects)
%           Usually, you will only use this option if you have run the
%           program before and wish to use a saved mask (so that you don't
%           have to repick the points for the regions of interest.
%               Subjects is a structure array containing individual Subject structures
%               with the following fields:
%                   Subject = Subjects{1};
%                   Subject.subject = {'Name'};
%                   Subject.allmasks = [];
%                   Subject.regions = [];  
%                   Subject.mask_perims = {};
%                   Subject.mask_perim_nums = [];
%
% 2. Modify parameters
%   a. timePoints for CDR and BL file
%   b. rotation angles- determined from Curry
%   c. cortex image name- a .jpeg file of a top-down view of the cortex
%   d. base directory
%   e. tasklist
%   f. name of .pom files to get points from Curry
%   g. 


function [OI,OI_interp,ActiveRatio,AOI_L,max_str,center,dis2censul]=target_ccd_area_BL2(subjects)

warning off MATLAB:divideByZero

Base_dir='C:\Albert Chen\Subject\';

%**************************************************************************

%1. ENTER SUBJECT PARAMETERS BEFORE RUNNING THE PROGRAM
%------------------------------------------------
% subject = {'CM1new'}
% 
% method={'Loreta1_0ICA'};
% %taskList={'abd','we','ws'};
% taskList={'abd'};
%     
% timePoints_total = [207, 104, 104];          %number of time points in each task
% timePoints_BL_total = [28, 52, 52];       %number of time points in each baseline
% timeRange_total = {[-700 100],[-300 100],[-300 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1950 -1850],[-1900 -1700],[-1900 -1700]};   %time ranges of baselines
% %%USING FAKE BASELINES
% %COORDINATES ARE ROTATED 21degrees about the x-axis
% alpha = -21*pi/180;
% beta = 0*pi/180;
% gamma = 0*pi/180;

% im_cortex_name = 'CMcortex2.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%------------------------------------------------
% subject = {'CM2new'}
% 
% method={'Loreta1_0ICA'};
% %taskList={'add','ef','hc','tab'};
% taskList={'ef'};
% 
% timePoints_total = [207, 104, 104, 104];          %number of time points in each task
% timePoints_BL_total = [52, 52, 52, 52];       %number of time points in each baseline
% timeRange_total = {[-700 100],[-300 100],[-300 100],[-300 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1900 -1700],[-1900 -1700],[-1900 -1700],[-1900 -1700]};   %time ranges of baselines
% %COORDINATES ARE ROTATED 21degrees about the x-axis
% alpha = -21*pi/180;
% beta = 0*pi/180;
% gamma = 0*pi/180;
% im_cortex_name = 'CMcortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%------------------------------------------------
% subject = {'CM042706'}
% 
% method={'Loreta1_0ICA'};
% %taskList={'abd_un','ef_un','hc_un'};   %list of tasks in this directory
% taskList={'abd_un'};
%     
% timePoints_total = [207, 104, 104];          %number of time points in each task
% timePoints_BL_total = [52, 52, 52];       %number of time points in each baseline
% timeRange_total = {[-700 100],[-300 100],[-300 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1900 -1700],[-1900 -1700],[-1900 -1700]};   %time ranges of baselines
% %COORDINATES ARE ROTATED 18degrees about the x-axis
% alpha = -21*pi/180;
% beta = 0*pi/180;
% gamma = 0*pi/180;
% im_cortex_name = 'CMcortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%------------------------------------------------
% subject = {'WSrobot'}
% subjects = subject;
% 
% %method={'Loreta1_0ICA','Loreta1_2ICA'};
% method={'Loreta1_5ICANEW'};
% %taskList={'RE,'zeroRE','twofiveRE'};
% taskList={'RE'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [28, 52];       %number of time points in each baseline
% timeRange_total = {[-700 100],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1950 -1800],[-1950 -1750]};   %time ranges of baselines
% 
% im_cortex_name = 'WSrobotcortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%-------------------------------------------------------------------
% subject = {'SNrobot'}
% 
% method={'Loreta1'};
% taskList={'25RE'};
%     
% timePoints_total = [207, 104, 104];          %number of time points in each task
% timePoints_BL_total = [41, 52, 52];       %number of time points in each baseline
% timeRange_total = {[-700 100],[-300 100],[-300 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1950 -1800],[-1900 -1700],[-1900 -1700]};   %time ranges of baselines
% %%USING FAKE BASELINES
% %COORDINATES ARE ROTATED 21degrees about the x-axis
% alpha = -21*pi/180;
% beta = 0*pi/180;
% gamma = 0*pi/180;

% im_cortex_name = 'SNcortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%------------------------------------------------

% subject = {'GMrobot'}
% 
% method={'Loreta1'};
% taskList={'25RE'};
%     
% timePoints_total = [207, 104, 104];          %number of time points in each task
% timePoints_BL_total = [3, 52, 52];       %number of time points in each baseline
% timeRange_total = {[-700 100],[-300 100],[-300 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1950 -1800],[-1900 -1700],[-1900 -1700]};   %time ranges of baselines
% %%USING FAKE BASELINES
% %COORDINATES ARE ROTATED 21degrees about the x-axis
% alpha = -21*pi/180;
% beta = 0*pi/180;
% gamma = 0*pi/180;

% im_cortex_name = 'GMcortex2.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT

%--------------------------------------------------------------------------
subject = {'AKpre1'}
subjects = subject;

%method={'Loreta1_0ICA','Loreta1_2ICA'};
method={''};
%taskList={'RE,'zeroRE','twofiveRE'};
taskList={'abd_un'};

timePoints_total = [207, 207];          %number of time points in each task
timePoints_BL_total = [41, 52];       %number of time points in each baseline
timeRange_total = {[-700 100],[-700 100]};     %time ranges of all tasks
timeRange_BL_total = {[-1950 -1800],[-1950 -1750]};   %time ranges of baselines

im_cortex_name = 'AKpre1cortex.jpg';
plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT


%*************************************************************************

%2. SOME INITIALIZATION

%individual time plots, what to plot in "for(timestart:timeinc:timestop)"
timestart=1;
timestop=max(timePoints_total);
timeinc=ceil((timestop-timestart+1)/16);  %want about 16 individual plots

%transformation matrix to coregister curry reconstruction file with cortex image
%COORDINATES ARE ROTATED 18degrees about the x-axis
alpha = -21*pi/180; 
beta = 0*pi/180;    
gamma = 0*pi/180; 

Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];    %rotation about x axis
Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];        %rotation about y axis
Rz = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];    %rotation about z axis
Rx = Rx*Ry*Rz;

%average cortex plot, what to average in plot in "for plotstart:plotstop"
%usually from -300ms to 0ms or something
%so have to figure out which indices correspond to those times
plotstart=1;

show_frame=0;
debug_yj=1;
subNo=length(subjects);
%**************************************************************************

%3. LOOP THROUGH ALL SUBJECTS
for k=1:subNo
    
    subjectName=subjects{k};    %subject name, name of directory where data is stored
    
    winPoints=15;
    
    %LOOP THROUGH ALL TASKS
    for i=1:length(taskList)
        cur_task=taskList{i};
        
        %----------------------------
        %modify number of timepoints used
        timePoints = timePoints_total(k);
        timePoints_BL = timePoints_BL_total(k);
        %----------------------------
        
        %LOOP THROUGH ALL METHODS
        for methodInd=1:length(method)
            ActCount=0;
%             cdr_file=[cur_task,'_',method{methodInd},'.cdr'];
%             cdr_BL_file=[cur_task,'_BL_',method{methodInd},'.cdr'];
            cdr_file=[cur_task,'.cdr'];
            cdr_BL_file=[cur_task,'_BL','.cdr'];
%             cdr_file=[cur_task,'.cdr'];
            
            disp('********************');
            disp(cdr_file)
            cdr_file_name=[Base_dir,subjectName,'\new\',cdr_file];
            cdr_BL_file_name=[Base_dir,subjectName,'\new\',cdr_BL_file];
            %--------------------------------------------------------------
            %read in the index of ROI and the rotation matrix
            if (i==1 & methodInd==1)
                
                %read in index, rotation matrix    
                [ROI_ind,Rx]=select_ROI_AC(subjectName,-21);
                
                %Read in central sulcus file, may be different name
                censul_file_name1=[Base_dir,subjectName,'\Points\','CS_cc.pom'];   %central sulcus
                
                %if only one central sulcus file
                [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name1,'LOCATION',0,0);
                censul2 = censul1;
                censul=[censul1];
                censul=(Rx*censul')';
            end  
            
            %--------------------------------------------------------------
            %get all cortex locations from *.cdr file
            fprintf(1, 'get cortex locations... ')
            [cortexL,Lcount,LNR]=read_Curry_file3(cdr_file_name,'LOCATION',0,0);
            fprintf(1, 'done\n')
            
            %get all cortex strengths corresponding closest to edge points
            fprintf(1, 'get cortex strengths... ')
            [cortexC,targetL,targetC]=find_cortex_strengths_Jun(cdr_file_name,cortexL,ROI_ind,timePoints,0);
            fprintf(1, 'done\n')
            
            %get baseline strengths of points in target areas
            [cortexC_BL,targetL_BL,targetC_BL]=find_cortex_strengths_Jun(cdr_BL_file_name,cortexL,ROI_ind,timePoints_BL,0);
            %--------------------------------------------------------------
            
            %rotate cortex and target points
            cortexL=(Rx*cortexL')';
            targetL=(Rx*targetL')';
            
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

            targetC_BL=mean(targetC_BL,2);
     
            %targetC takes last number of points equal to winPoints
%             targetC=targetC(:,(end-winPoints+1):end);
            winPoints = plotstop-plotstart+1;
            targetC=targetC(:,plotstart:plotstop);            

            %subtract baseline
%             targetC=targetC-repmat(targetC_BL,1,winPoints);
            targetC=targetC-repmat(targetC_BL,1,plotstop-plotstart+1);
            targetC=targetC(:,plotstart:plotstop);

            [locationNo,TimeNo]=size(targetC);
            
            targetC_sum=sum(targetC,2);
            targetC_norm=targetC_sum/max(targetC_sum);
            
            
            %--------Jun added on 09/25/06 for the case that some
            %targetC_BL is larger than the targetC
            ind=find(targetC_norm<0);
            targetC_norm(ind)=0;
            %--------------------------------------------------------
            
            %-------Jun added for looking for the CoG---------------
            center(k,i,:)=sum(repmat(targetC_norm,1,3).*targetL,1);
            center(k,i,:)=center(k,i,:)/sum(targetC_norm);
            max_str(k,i)=max(targetC_norm);
            [dis2censul(k,i),index]=min(sqrt( sum( (censul-repmat(squeeze(center(k,i,:))',size(censul,1),1)).^2,2 )  ));

            if debug_yj
                figure
                hold on
                CI = griddata(targetL(:,1),targetL(:,2),targetC_sum/winPoints,XI,YI);
                surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
                h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
                h2=plot3(center(k,i,1),center(k,i,2),center(k,i,3),'m+');
                view(90,90)
                colorbar
                
%                 pause
            end
            
            %--------------------------------------------------------------
            %FIND ACTIVE RATIO FOR ENTIRE MOTOR AREA
            
%             phat_exp(k,i)=mle(targetC_norm,'distribution', 'Exponential')
%             phat_norm(k,i)=mle(targetC_norm,'distribution','Normal')
            
            [muhat, muci] = expfit(targetC_norm);
            ind=find(targetC_norm>muci(2));
            
            ActiveRatio(k,i)=length(ind)/locationNo; %active area ratio
            
            %--------------------------------------------------
            %compare to Albert's method of finding active ratio
            %--------------------------------------------------------------
            
            
            
        end %end of methodInd
        CCD_task(:,i)=targetC_sum;
        CCD_task_norm(:,i)=targetC_norm;
    end % end of task
    
    
    OverArea_norm=CCD_task_norm(:,1).*CCD_task_norm(:,2);
    OverArea=CCD_task(:,1).*CCD_task(:,2);
    OI(k,1)=sum(OverArea)/sum(CCD_task(:,2)); %overlap index counting the strength
    OI(k,2)=sum(OverArea)/sum(CCD_task(:,1));
    CI_task1 = griddata(targetL(:,1),targetL(:,2),CCD_task_norm(:,1),XI,YI);
    CI_task2 = griddata(targetL(:,1),targetL(:,2),CCD_task_norm(:,2),XI,YI);
    TA=isnan(CI_task1);
	ind3=find(TA==1);
	CI_task1(ind3)=0;
    TA=isnan(CI_task2);
	ind3=find(TA==1);
	CI_task2(ind3)=0;
   
    OverArea_inter= CI_task1.* CI_task2;

    OI_interp(k,1)=sum(sum(OverArea_inter,1),2)/sum(sum(CI_task2,1),2); %overlap index counting the strength and interplate
    OI_interp(k,2)=sum(sum(OverArea_inter,1),2)/sum(sum(CI_task1,1),2);

    [muhat, muci] = expfit(OverArea_norm);
    ind_overlap=find(OverArea_norm>muci(2));
    [muhat, muci] = expfit(CCD_task_norm(:,1));
    ind1=find(CCD_task_norm(:,1)>muci(2));
    [muhat, muci] = expfit(CCD_task_norm(:,2));
    ind2=find(CCD_task_norm(:,2)>muci(2));
    AOI_L(k)=2*length(ind_overlap)/(length(ind1)+length(ind2)); %modified overlap area index
%     OverArea_s(k)=sum(OverArea_norm)/locationNo;    

%     [muhat,sigmahat,muci,sigmaci] = normfit(OverArea_norm);
%     ind_overlap=find(OverArea_norm>muci(2));
%     [muhat,sigmahat,muci,sigmaci] = normfit(CCD_task_norm(:,1));
%     ind1=find(CCD_task_norm(:,1)>muci(2));
%     [muhat,sigmahat,muci,sigmaci] = normfit(CCD_task_norm(:,2));
%     ind2=find(CCD_task_norm(:,2)>muci(2));
%     AOI_L(k)=2*length(ind_overlap)/(length(ind1)+length(ind2)); %modified overlap area index


    
    if show_frame
        figure 
        CI = griddata(targetL(:,1),targetL(:,2),OverArea,XI,YI);
        hold on
        surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)
        h=plot3(censul(:,1),censul(:,2),censul(:,3),'g.');
        h2=plot3(center(:,1),center(:,2),center(:,3),'g.');
        colorbar
        set(h,'LineWidth',7)
        view(90,90)
    end

    CI = griddata(targetL(:,1),targetL(:,2),OverArea_norm,XI,YI);
	TA=isnan(CI);
	ind3=find(TA==1);
	CI(ind3)=0;

	[counts,xx]=imhist(CI);
    [muhat, muci] = expfit(counts);
    ind=find(counts>muci(2));
    thr=xx(ind(end));
    ind2=find(CI>thr);
    AOI_CI(k)=length(ind2);
    ind4=find(CI>0);
    LNo_CI(k)=length(ind4);
    AOI_inter(k)=AOI_CI(k)/LNo_CI(k);

    clear CCD_task CCD_task_norm OverArea OverArea_norm XI YI ZI CI
end % end of subject Number

beep


