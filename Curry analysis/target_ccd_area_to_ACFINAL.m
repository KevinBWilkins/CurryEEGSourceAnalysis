% TARGET_CCD_AREA_TO_AC()
% function ActiveRatio=target_ccd_area_to_ACFINAL(subject)
% By Albert Chen, Jun Yao
%
% This function read the cortical current density (CCD) information
% in the target local areas and then finds the size of the active area
% in each region of interest. It can use either ANOVA or arbitrary 
% thresholding to determine the active area.
%
% Last Modified 9/13/06 - AC
%
% Manual changes must be made in the following sections:
%   1. Enter subject parameters before running program, specify when 0ms occurs
%
%
% Sample usage:
% 
% target_ccd_area_to_ACFINAL()
% [Subject] = target_ccd_area_to_ACFINAL(Subject)


function [Subject] = target_ccd_area_to_ACFINAL(varargin)

if length(varargin) > 0
    %load previous values
    Subject = varargin{1};
    maxvalues = Subject.maxvalues;
    newThresholds = Subject.newThresholds;
    percentmaxthresholds = Subject.percentmaxthresholds;
    SIs = Subject.SIs;
    AAIs = Subject.AAIs;
    allmasks = Subject.allmasks;
    if isfield(Subject,'mask_perims')
        mask_perims = Subject.mask_perims;
        mask_perim_nums = Subject.mask_perim_nums;
        regions = Subject.regions;
    else
        mask_perims = {};
        mask_perim_nums = [];
    end
else
    mask_perim_nums = [];
end


Base_dir='C:\Albert Chen\Subject\';
warning off MATLAB:divideByZero

%1. Enter subject parameters before running the program

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
% Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
% Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
% Rz = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
% 
% Rx = Rx*Ry*Rz;
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
% Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
% Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
% 
% Rx = Rx*Ry;
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
% Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
% Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
% 
% Rx = Rx*Ry;
% im_cortex_name = 'CMcortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%------------------------------------------------
% subject = {'BHrobot'}
% 
% method={'Loreta1_2ICANEW'}
% taskList={'RE2'}
% 
% timePoints_total = [207, 90, 90];          %number of time points in each task
% timePoints_BL_total = [52, 26, 26];       %number of time points in each baseline
% timeRange_total = {[-700 50],[-300 50],[-300 50]};     %time ranges of all tasks
% timeRange_BL_total = {[-1950 -1800],[-1900 -1800],[-1900 -1800]};   %time ranges of baselines
% %COORDINATES ARE ROTATED 18degrees about the x-axis
% alpha = -35*pi/180;
% beta = 7*pi/180;
% Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
% Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
% 
% Rx = Rx*Ry;
% im_cortex_name = 'BHcortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%------------------------------------------------
% subject = {'BHrobot'}
% 
% method={'Loreta1'}
% taskList={'RE'}
% 
% timePoints_total = [90, 90, 90];          %number of time points in each task
% timePoints_BL_total = [52, 26, 26];       %number of time points in each baseline
% timeRange_total = {[-300 50],[-300 50],[-300 50]};     %time ranges of all tasks
% timeRange_BL_total = {[-1900 -1700],[-1900 -1800],[-1900 -1800]};   %time ranges of baselines
% %COORDINATES ARE ROTATED 35degrees about the x-axis, 7degrees around y-axis
% alpha = -35*pi/180;
% beta = 7*pi/180;
% Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
% Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
% 
% Rx = Rx*Ry;
% im_cortex_name = 'BHcortex.jpg';
% plotstop=78;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%------------------------------------------------
% subject = {'PGrobot'}
% 
% method={'Loreta1'};
% % taskList={'newREtwofive','newREzero','newRE'};
% taskList = {'newREtwofive'};
% 
% timePoints_total = [91, 91, 91];          %number of time points in each task
% timePoints_BL_total = [52, 52, 52];       %number of time points in each baseline
% timeRange_total = {[-300 50],[-300 50],[-300 50]};     %time ranges of all tasks
% timeRange_BL_total = {[-1900 -1700],[-1900 -1700],[-1900 -1700]};   %time ranges of baselines

%------------------------------------------------
subject = {'WSrobot'}

%method={'Loreta1_0ICA','Loreta1_2ICA'};
method={'Loreta1_5ICANEW'};
%taskList={'RE,'zeroRE','twofiveRE'};
taskList={'RE'};

timePoints_total = [207, 207];          %number of time points in each task
timePoints_BL_total = [28, 52];       %number of time points in each baseline
timeRange_total = {[-700 100],[-700 100]};     %time ranges of all tasks
timeRange_BL_total = {[-1950 -1850],[-1950 -1750]};   %time ranges of baselines

%COORDINATES ARE ROTATED 18degrees about the x-axis
alpha = -21*pi/180;
beta = 0*pi/180;
gamma = 0*pi/180;
Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
Rz = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];

Rx = Rx*Ry*Rz;
im_cortex_name = 'WScortex2.jpg';
plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%-------------------------------------------------------------------
% subject = {'VOrobot2'}
% 
% method={'Loreta1'}
% taskList={'RE'}
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-700 100],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1950 -1800],[-1950 -1750]};   %time ranges of baselines
% %COORDINATES ARE ROTATED 18degrees about the x-axis
% alpha = -21*pi/180;
% beta = 0*pi/180;
% Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
% Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
% 
% Rx = Rx*Ry;
% im_cortex_name = 'VOcortex2.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT

%-------------------------------------------------------------------
% subject = {'BHrobot2'}
% 
% method={'Loreta1'}
% taskList={'RE'}
% 
% timePoints_total = [155, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-700 -100],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1950 -1800],[-1950 -1750]};   %time ranges of baselines
% %COORDINATES ARE ROTATED 18degrees about the x-axis
% alpha = -3*pi/180;      %positive goes up, negative goes down
% beta = -2*pi/180;       %positive goes left, negative goes right
% Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
% Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
% 
% Rx = Rx*Ry;
% im_cortex_name = 'BHcortex2.jpg';
% plotstop=155;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%---------------------------------------------------------------------
% 
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
% Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
% Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
% Rz = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
% 
% Rx = Rx*Ry*Rz;
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
% Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
% Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
% Rz = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
% 
% Rx = Rx*Ry*Rz;
% im_cortex_name = 'GMcortex2.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT


%*************************************************************************
%2. Some initialization

%individual time plots, what to plot in "for(timestart:timeinc:timestop)"
timestart=1;
timestop=max(timePoints_total);
timeinc=ceil((timestop-timestart+1)/16);  %want about 16 individual plots

%average cortex plot, what to average in plot in "for plotstart:plotstop"
%usually from -300ms to 0ms or something
%so have to figure out which indices correspond to those times
plotstart=1;

%**************************************************************************



subNo=length(subject);
phaseList={'whole'};
plotCort=0;
show_frame=0;
yjcolor=['b','r','m'];
currents={};
V1='B'; %baseline
V2='E'; %activation (execution)


%3. Run through loop of subjects

for k=1:subNo
    subjectName=subject{k};
    
    
    %--------------------------------------------------------------------
    %4. Import areas of interest
    
    %using Curry 4.0
%     edge1_file_name=[Base_dir,subjectName,'\Data\MRI\','0002.sp8'];M1_rt
%     edge2_file_name=[Base_dir,subjectName,'\Data\MRI\','0002.sp9'];M1_lt
%     edge3_file_name=[Base_dir,subjectName,'\Data\MRI\','0002.sp6'];S1_rt
%     edge4_file_name=[Base_dir,subjectName,'\Data\MRI\','0002.sp7'];S1_lt
%     edge5_file_name=[Base_dir,subjectName,'\Data\MRI\','0002.sp2'];PM_rt
%     edge6_file_name=[Base_dir,subjectName,'\Data\MRI\','0002.sp3'];PM_lt
%     edge7_file_name=[Base_dir,subjectName,'\Data\MRI\','0002.sp0'];SMA_rt
%     edge8_file_name=[Base_dir,subjectName,'\Data\MRI\','0002.sp1'];SMA_lt
	
    %using Curry 5.0
%     edge1_file_name=[Base_dir,subjectName,'\Data\MRI\','M1_rt.pom'];%M1_rt
%     edge2_file_name=[Base_dir,subjectName,'\Data\MRI\','M1_lt.pom'];%M1_lt
%     edge3_file_name=[Base_dir,subjectName,'\Data\MRI\','S1_rt.pom'];%S1_rt
%     edge4_file_name=[Base_dir,subjectName,'\Data\MRI\','S1_lt.pom'];%S1_lt
% 	edge5_file_name=[Base_dir,subjectName,'\Data\MRI\','PM_rt.pom'];%PM_rt
%     edge6_file_name=[Base_dir,subjectName,'\Data\MRI\','PM_lt.pom'];%PM_lt
%     edge7_file_name=[Base_dir,subjectName,'\Data\MRI\','SMA_rt.pom'];%SMA_rt
%     edge8_file_name=[Base_dir,subjectName,'\Data\MRI\','SMA_lt.pom'];%SMA_lt

    %using Curry 5.0 Carolina's pts
    edge1_file_name=[Base_dir,subjectName,'\Data\MRI\','M1_rt_cc.pom'];%M1_rt
    edge2_file_name=[Base_dir,subjectName,'\Data\MRI\','M1_lt_cc.pom'];%M1_lt
    edge3_file_name=[Base_dir,subjectName,'\Data\MRI\','S1_rt_cc.pom'];%S1_rt
    edge4_file_name=[Base_dir,subjectName,'\Data\MRI\','S1_lt_cc.pom'];%S1_lt
	edge5_file_name=[Base_dir,subjectName,'\Data\MRI\','PM_rt_cc.pom'];%PM_rt
    edge6_file_name=[Base_dir,subjectName,'\Data\MRI\','PM_lt_cc.pom'];%PM_lt
    edge7_file_name=[Base_dir,subjectName,'\Data\MRI\','SMA_rt_cc.pom'];%SMA_rt
    edge8_file_name=[Base_dir,subjectName,'\Data\MRI\','SMA_lt_cc.pom'];%SMA_lt

    
    %READ IN LOCATIONS OF POINTS
    [edge1,Ecount,ENR]=read_Curry_file4_AC(edge1_file_name,'LOCATION',0,0);
	[edge2,Ecount,ENR]=read_Curry_file4_AC(edge2_file_name,'LOCATION',0,0);
	[edge3,Ecount,ENR]=read_Curry_file4_AC(edge3_file_name,'LOCATION',0,0);
	[edge4,Ecount,ENR]=read_Curry_file4_AC(edge4_file_name,'LOCATION',0,0);
    [edge5,Ecount,ENR]=read_Curry_file4_AC(edge5_file_name,'LOCATION',0,0);
	[edge6,Ecount,ENR]=read_Curry_file4_AC(edge6_file_name,'LOCATION',0,0);
	[edge7,Ecount,ENR]=read_Curry_file4_AC(edge7_file_name,'LOCATION',0,0);
	[edge8,Ecount,ENR]=read_Curry_file4_AC(edge8_file_name,'LOCATION',0,0);
%     [edge9,Ecount,ENR]=read_Curry_file4_AC(edge9_file_name,'LOCATION',0,0);
%     [edge10,Ecount,ENR]=read_Curry_file4_AC(edge10_file_name,'LOCATION',0,0);

       
    %--------------------------------------------
    %5. Read in central sulcus file, may be different name
    % also read in diamond or pentagon for reference if necessary
    
    %curry 4.0
%     censul_file_name1=[Base_dir,subjectName,'\Data\MRI\','0002.sp4'];   %central sulcus_rt
%     censul_file_name2=[Base_dir,subjectName,'\Data\MRI\','0002.sp5'];   %central sulcus_lt
     
    %curry 5.0
%      censul_file_name1=[Base_dir,subjectName,'\Data\MRI\','0002.sp0'];   %central sulcus_rt
%      censul_file_name2=[Base_dir,subjectName,'\Data\MRI\','0002.sp1'];   %central sulcus_lt

%     censul_file_name1=[Base_dir,subjectName,'\Data\MRI\','centralsulcus_rt.pom'];   %central sulcus_rt
%     censul_file_name2=[Base_dir,subjectName,'\Data\MRI\','centralsulcus_lt.pom'];   %central sulcus_lt

    censul_file_name1=[Base_dir,subjectName,'\Data\MRI\','CS_cc.pom'];   %central sulcus
    
    %if two central sulcus files
%     [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name1,'LOCATION',0,0);
%     [censul2,Ecount,ENR]=read_Curry_file3(censul_file_name2,'LOCATION',0,0);
%     censul=[censul1;censul2];
    
    %if only one central sulcus file
    [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name1,'LOCATION',0,0);
    censul2 = censul1;
    censul=[censul1];

%     diam_file_name1=[Base_dir,subjectName,'\Data\MRI\','0003.sp4'];   %diamond or pentagon for reference
%     [diam1,Ecount,ENR]=read_Curry_file3(diam_file_name1,'LOCATION',0,0);    
    
    %----------------------------------------------
    %6. Combine right and left sides of areas to make whole SMA, PM, M1, and S1
%     areanames = {'SMA','PM','M1','S1','ContraHem','IpsiHem','Whole'};
    areanames = {'SMA_rt','PM_rt','M1_rt','S1_rt','ContraHem','IpsiHem','Whole'};
%     areanames = {'SMA_lt','PM_lt','M1_lt','S1_lt','ContraHem','IpsiHem','Whole'};
    
    areas{1} = [edge7; edge8]; %SMA
    areas{2} = [edge5; edge6]; %PM
    
    areas{3} = [edge1; edge2]; %M1
    areas{4} = [edge3; edge4]; %S1
    
%     areas{1} = [edge8];   %Ipsilateral areas only
%     areas{2} = [edge6];
%     areas{3} = [edge2];
%     areas{4} = [edge4];
    
    areas{1} = [edge7];     %Contralateral areas only
    areas{2} = [edge5];
    areas{3} = [edge1];
    areas{4} = [edge3];
    
    areas{5} = [edge1;edge3;edge5;edge7];       %RIGHT HEMISPHERE
    areas{6} = [edge2;edge4;edge6;edge8];       %LEFT HEMISPHERE
%     areas{7} = [edge1;edge3;edge5;edge7;edge2;edge4;edge6;edge8];   %ALL AREAS
    areas{7} = [edge1;edge5;edge7;edge2;edge6;edge8;edge3;edge4];   %ALL MOTOR AREAS
    
%     areas{7} = [edge1;edge3;edge5;edge7;edge2;edge4;edge6;edge8;edge9;edge10];   %ALL AREAS including PPC
    
    %end of picking areas of interest
    %---------------------------------------------------------------
    
       
    %----------------------------------------------
    %7. Loop through tasks- such as RE, 0RE, 25RE, etc
    for i=1:length(taskList)
        cur_task=taskList{i};
        
        timePoints = timePoints_total(i);
        timePoints_BL = timePoints_BL_total(i);
        timeRange = timeRange_total{i};
        timeRange_BL = timeRange_BL_total{i};
        
        %---------------------------------------------------------------
        %8. Loop through areas of interest- 1=M1, 2=S1, 3=PM, 4=SMA,
        %5=right hem, 6=left hem, 7=whole cortex
        for areanum = 7:7
            fprintf(1, 'processing areanum %d... \n',areanum)
            edge = areas{areanum};
            
            %----------------------------------
            %9. Loop through methods- Loreta1, etc
            for methodInd=1:length(method)
                ActCount=0;
                
                cdr_file=[cur_task,'_',method{methodInd},'.cdr'];
                cdr_BL_file=[cur_task,'_baseline_',method{methodInd},'.cdr'];
                
                disp('********************');
                disp(cdr_file)
                cdr_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_file];
                cdr_BL_file_name=[Base_dir,subjectName,'\results_inv\noVMOV\',cdr_BL_file];%BASELINE
                
                %--------------------------------------------
                %10. Get strengths of points in target areas
                %use 1 for last value to display picture
                
                %get all cortex locations
                fprintf(1, 'get cortex locations... ')
                [cortexL,Lcount,LNR]=read_Curry_file3(cdr_file_name,'LOCATION',0,0);
                fprintf(1, 'done\n')
                
                %get all cortex strengths corresponding closest to edge points
                fprintf(1, 'get cortex strengths... ')
                [cortexC,targetL,targetC]=find_cortex_strengths_AC(cdr_file_name,cortexL,edge,timePoints,0);
                %[cortexL,targetL,targetC]=find_target_PC3_AC(cdr_file_name,cortexL,edge,timePoints,0);
                fprintf(1, 'done\n')
                
                %11. Get baseline strengths of points in target areas
                [cortexC_BL,targetL_BL,targetC_BL]=find_cortex_strengths_AC(cdr_BL_file_name,cortexL,edge,timePoints_BL,0);
                
                %take first 50 ms of timerange and use as baseline if
                %didn't run separate baseline
                totaltime = timeRange(2)-timeRange(1);
                times1 = timeRange(1):totaltime/(timePoints-1):timeRange(1)+50;
                totaltime_BL = timeRange_BL(2)-timeRange_BL(1);
                times1_BL = timeRange_BL(1):totaltime_BL/(timePoints_BL-1):timeRange_BL(2);
                times1_seg = timeRange(1):totaltime/(timePoints-1):timeRange(1)+50;
                meantargetC_BL = mean(targetC_BL,1);
                totalmeantargetC_BL = mean(meantargetC_BL);
                
                
                [m,n] = size(targetC_BL);
                longtargetC_BL = reshape(targetC_BL,m*n,1);
                sorttargetC_BL = sort(longtargetC_BL);
                stdtargetC_BL = std(sorttargetC_BL);
                
                [m,n] = size(targetC);
                longtargetC = reshape(targetC,m*n,1);
                sorttargetC = sort(longtargetC);
                
                [m,n] = size(targetC(:,1:length(times1_seg)));
                longtargetC_seg = reshape(targetC(:,1:length(times1_seg)),m*n,1);
                sorttargetC_seg = sort(longtargetC_seg);
                
                maxtargetC = max(max(targetC(:,1:plotstop)));
                                
                %thresholds are a percentage of maximum value inside that
                %area   %DOES NOT USE SINGLE SET OF DATA FOR MAX VALUE TO STANDARDIZE
                %ACROSS DIFFERENT SETS OF DATA FOR THE SAME PERSON

                %Threshold = totalmeantargetC_BL + 2*stdtargetC_BL;
                Threshold = max(sorttargetC_BL)
                
                thresholds = [.85; .80; .75; .70; .65; .60; .50; .40; .20]*maxtargetC;
                
                PercentThresholdofMax = Threshold/maxtargetC
                
                %Plot of thresholds
%                 figure(80+i)
%                 plot(sorttargetC_BL,'Color','g','Linewidth',2)
%                 hold on;
%                 plot(sorttargetC,'Color','b','Linewidth',2)
%                 plot(sorttargetC_seg,'Color','k','Linewidth',2)
%                 %line([1 length(sorttargetC)],[totalmeantargetC_BL totalmeantargetC_BL],'Color','r','Linewidth',2)
%                 %line([1 length(sorttargetC)],[totalmeantargetC_BL+2*stdtargetC_BL totalmeantargetC_BL+2*stdtargetC_BL],'Color','r','Linewidth',2)
%                 %line([1 length(sorttargetC)],[totalmeantargetC_BL-2*stdtargetC_BL totalmeantargetC_BL-2*stdtargetC_BL],'Color','r','Linewidth',2)
%                 line([1 length(sorttargetC)],[Threshold Threshold],'Color','r','Linewidth',2)
%                 line([1 length(sorttargetC)],[thresholds thresholds],'Color','m','Linewidth',2)
                
                thresholds = Threshold;
                
                
                
                %11b. Subtract baseline strengths from all strengths
                %targetC = targetC - repmat(meantargetC_BL,1,timePoints);
                
                %------------------------------------------
                %12. Plot pretty pictures
                samp_inc = 0.8;
                
                %target locations- unrotated
                minX=min(targetL(:,1));
                maxX=max(targetL(:,1));
                minY=min(targetL(:,2));
                maxY=max(targetL(:,2));
                minZ=min(targetL(:,3));
                maxZ=max(targetL(:,3));
                x=[minX:samp_inc:maxX];
                y=[minY:samp_inc:maxY];
                [XI,YI] = meshgrid(x,y); 
                ZI = griddata(targetL(:,1),targetL(:,2),targetL(:,3),XI,YI);
                
                %ROTATE ALL TARGET LOCATIONS before create mask
                %find limits of rotated target locations
                rot_targetL = (Rx*targetL')';
                minrot_X=min(rot_targetL(:,1));
                maxrot_X=max(rot_targetL(:,1));
                minrot_Y=min(rot_targetL(:,2));
                maxrot_Y=max(rot_targetL(:,2));
                minrot_Z=min(rot_targetL(:,3));
                maxrot_Z=max(rot_targetL(:,3));
                rot_x=[minrot_X:samp_inc:maxrot_X];
                rot_y=[minrot_Y:samp_inc:maxrot_Y];
                [rot_XI,rot_YI] = meshgrid(rot_x,rot_y); 
                rot_ZI = griddata(rot_targetL(:,1),rot_targetL(:,2),rot_targetL(:,3),rot_XI,rot_YI);
                
                %ROTATE ALL EDGE LOCATIONS before create mask
                %find limits of rotated edge locations
                rot_edgeL = (Rx*edge')';
                minrot_edgeX=min(rot_edgeL(:,1));
                maxrot_edgeX=max(rot_edgeL(:,1));
                minrot_edgeY=min(rot_edgeL(:,2));
                maxrot_edgeY=max(rot_edgeL(:,2));
                minrot_edgeZ=min(rot_edgeL(:,3));
                maxrot_edgeZ=max(rot_edgeL(:,3));
                rot_edgex=[minrot_edgeX:samp_inc:maxrot_edgeX];
                rot_edgey=[minrot_edgeY:samp_inc:maxrot_edgeY];
                [rot_edgeXI,rot_edgeYI] = meshgrid(rot_edgex,rot_edgey); 
                rot_edgeZI = griddata(rot_edgeL(:,1),rot_edgeL(:,2),rot_edgeL(:,3),rot_edgeXI,rot_edgeYI);
                
                %ROTATE CORTEX LOCATIONS TO MATCH TOP DOWN LOOK
                rot_cortexL = (Rx*cortexL')';
                mincortX=min(rot_cortexL(:,1));
                maxcortX=max(rot_cortexL(:,1));
                mincortY=min(rot_cortexL(:,2));
                maxcortY=max(rot_cortexL(:,2));
                bigx=[mincortX:samp_inc:maxcortX];
                bigy=[mincortY:samp_inc:maxcortY];
                l_bigy = length(bigy);      %dimensions of big cortex locations
                l_bigx = length(bigx);
                
                [bigXI,bigYI] = meshgrid(bigx,bigy);
                
                %get highest z-values at each x-y coordinate (only top surface of cortex)
                
                indzs = [];
                maxzs = [];
                %initial bigZI
                bigZI = griddata(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),bigXI,bigYI);
                for ind_x = 1:length(bigx)
                    [ind_x_found] = find((abs(rot_cortexL(:,1)-bigx(ind_x)))<samp_inc);
                    for ind_y = 1:length(bigy)
                        [ind_y_found] = find(abs((rot_cortexL(:,2)-bigy(ind_y)))<samp_inc);
                        if ~isempty(ind_y_found)
                            intersect_x_y = intersect(ind_x_found,ind_y_found);
                            if ~isempty(intersect_x_y)
                                
                                [maxz,indz] = max(rot_cortexL(intersect_x_y,3));
                                if bigx(ind_x)>=minrot_edgeX && bigx(ind_x)<=maxrot_edgeX && bigy(ind_y)>=minrot_edgeY && bigy(ind_y)<=maxrot_edgeY
                                    x_ind = ceil((bigx(ind_x)-minrot_edgeX)/samp_inc);
                                    y_ind = ceil((bigy(ind_y)-minrot_edgeY)/samp_inc);
                                    
                                    if y_ind<=length(rot_ZI(:,1)) && x_ind<=length(rot_ZI(1,:))
                                        newmaxz = rot_ZI(y_ind,x_ind);
                                    else
                                        newmaxz = NaN;
                                    end
                                    if ~isnan(newmaxz)
                                        indzs = [indzs; intersect_x_y(indz)];
                                        maxzs = [maxzs; newmaxz];
                                        bigZI(ind_y,ind_x) = maxzs(end);
                                    else
                                        indzs = [indzs; intersect_x_y(indz)];
                                        maxzs = [maxzs; maxz];
                                        bigZI(ind_y,ind_x) = maxzs(end);
                                    end
                                else
                                    indzs = [indzs; intersect_x_y(indz)];
                                    maxzs = [maxzs; maxz];
                                    bigZI(ind_y,ind_x) = maxzs(end);
                                end
                            end
                        end
                    end
                end
                [indzs,i_indzs,j_indzs] = unique(indzs);
                maxzs = maxzs(i_indzs);
                top_rot_cortexL = [rot_cortexL(indzs,1:2) maxzs];
                top_cortexC = cortexC(indzs,:);
                top_cortexC_BL = cortexC_BL(indzs,:);
                                
                
%                 figure
%                 plot3(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),'LineWidth',7,'color','r','LineStyle','.')
%                 hold on;
%                 plot3(top_rot_cortexL(:,1),top_rot_cortexL(:,2),top_rot_cortexL(:,3),'LineWidth',7,'color','b','LineStyle','.')
%                 plot3(rot_targetL(:,1),rot_targetL(:,2),rot_targetL(:,3),'LineWidth',7,'color','g','LineStyle','.')
%                 view(0,0);
                %bigZI = griddata(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),bigXI,bigYI);
                %bigZI = griddata(top_rot_cortexL(:,1),top_rot_cortexL(:,2),top_rot_cortexL(:,3),bigXI,bigYI);
                
                %ROTATE everything else too
                rot_censul = (Rx*censul')';
                  
                %------------------------------------------
                %12b. Find grid locations for sampled points on total grid
                %and create a mask where 1 = region, 0 = outside region
                
                %determine whether mask has already been created
                if length(varargin) < 1
                    createmaskflag = 1;
                else
                    size_allmasks = size(allmasks);
                    if size_allmasks(2) < areanum
                        createmaskflag = 1;
                    else
                        if isempty(allmasks{i,areanum,methodInd})
                            createmaskflag = 1;
                        else
                            createmaskflag = 0;
                        end
                    end
                end
                
                %IMAGE PROCESSING OF CORTEX PICTURE
                %used to overlay region of interest onto picture of cortex
                im_cortex_file_name=[Base_dir,subjectName,'\Data\MRI\',im_cortex_name];%image file location
                imA = imread(im_cortex_file_name);       %raw picture
                imB = im2bw(imA,.95);               %make black and white
                se = strel('square',5);             
                imC = imclose(imB,se);              %close up holes
                imD = ~imC;                         %flip white and black
                imE = imfill(imD,'holes');          %fill holes
                [im_i,im_j] = find(imE>0);          %find filled area
                imF = imE(min(im_i):max(im_i),min(im_j):max(im_j));     %
                imG = bwperim(imF,8);
                
                imA_resize = imA(min(im_i):max(im_i),min(im_j):max(im_j),:);
                %fit cortex picture to cortex pt locations
                
                [cortex_m,cortex_n] = size(imG);
                imA_resize2 = imresize(imA_resize,[cortex_m,cortex_n],'bilinear');
                
                %size of pixels in image
                samp_incX = (maxcortX-mincortX)/(cortex_n-1);
                samp_incY = (maxcortY-mincortY)/(cortex_m-1);
                
                %all x and y locations on the image
                x_locs = mincortX:samp_incX:maxcortX;
                y_locs = maxcortY:-samp_incY:mincortY;
                
                %find perimeter points (plotted in magenta)
                [imG_pts_i,imG_pts_j] = find(imG>0);
                for ind_cort=1:5:length(imG_pts_i)
                    loc_cortex(ind_cort,:) = [x_locs(imG_pts_j(ind_cort)),y_locs(imG_pts_i(ind_cort))];
                end
                loc_cortex = loc_cortex(1:5:end,:);                
                
                %determine midline
                %midline1_x = x_locs(round(length(x_locs)/2));
                [P,S] = polyfit(loc_cortex(:,2),loc_cortex(:,1),1);
                midline1_x = P(1)*y_locs+P(2);
                %[P,S] = polyfit(rot_cortexL(:,2),rot_cortexL(:,1),1)
                %midline2_x = P(1)*y_locs+P(2);
                
                %Figure displaying cortex picture, border of cortex, central sulcus, midline                
                figure
                imshow([mincortX maxcortX],[maxcortY mincortY],imA_resize2)
                axis on
                axis([mincortX maxcortX mincortY maxcortY])
                axis xy
                hold on;
                scatter3(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),7,'y');
                scatter3(loc_cortex(:,1),loc_cortex(:,2),zeros(length(loc_cortex(:,1)),1),7,'m')
                plot3(rot_censul(:,1),rot_censul(:,2),rot_censul(:,3),'Markersize',15,'color','b','LineStyle','.')
                %scatter3(rot_censul(:,1),rot_censul(:,2),rot_censul(:,3),3,'g');
                %scatter3(rot_targetL(:,1),rot_targetL(:,2),rot_targetL(:,3),7,'r');
                
                %plot midline
                %line([midline1_x midline1_x],[mincortY maxcortY],'Color','red','Linestyle','--')
                plot(midline1_x,y_locs,'color','m','Linewidth',2,'Linestyle','--')
                %plot(midline2_x,y_locs,'color','m','Linewidth',2,'Linestyle','--')
                
                view(0,90)
                axis on
                axis equal
                axis([mincortX maxcortX mincortY maxcortY])
                
                %resize imA_resize2 so that pixel size = samp_inc
                imA_resize3 = imresize(imA_resize2,[length(bigy),length(bigx)],'bilinear');
                
                
                %CREATE MASK OF ROTATED TARGET LOCATIONS
                if (createmaskflag)
                    fprintf(1, 'create mask... ')
                    
                    %create mask array- limits of ROTATED target locations                    
                    mask = zeros(length(rot_y),length(rot_x));
                    [mask_m,mask_n] = size(mask);
                    edge_mask = zeros(mask_m,mask_n);
                    
                    Xmat = repmat(rot_x,mask_m,1);
                    Ymat = repmat(rot_y,mask_n,1)';
                    
                    for loc_ind=1:length(rot_targetL(:,1))
                        %dx = rot_targetL(loc_ind,1); 
                        %dy = rot_targetL(loc_ind,2);
                        [minval,minrow] = min((rot_targetL(loc_ind,1)-Xmat).^2 + (rot_targetL(loc_ind,2)-Ymat).^2);
                        [minval,mincol] = min(minval);
                        
                        loc_y = minrow(mincol); 
                        loc_x = mincol;
                        %dx2 = rot_x(loc_x); 
                        %dy2 = rot_y(loc_y);
                        
                        mask(loc_y,loc_x) = 1;
                    end
                    
                    for loc_ind=1:length(rot_edgeL(:,1))
                        [minval,minrow] = min((rot_edgeL(loc_ind,1)-Xmat).^2 + (rot_edgeL(loc_ind,2)-Ymat).^2);
                        [minval,mincol] = min(minval);
                        
                        loc_y = minrow(mincol); 
                        loc_x = mincol;
                        edge_mask(loc_y,loc_x) = 1;
                    end
                    
                    mask2 = imclose(mask,se);
                    mask3 = imfill(mask2,'holes');
                    mask4 = imopen(mask3,se);
                    mask5 = mask4+2*edge_mask;   %show proposed pts with original edge pts selected in Curry
                    
                    
                    %DRAW CORTEX in background of mask selection
                    figure(100)
                    h1 = imagesc(imA_resize3);
                    %set(h1,'AlphaData',1);
                    
                    %make new mask with same size as imA_resize3
                    %figure out locations associated with imA_resize3
                    %figure out where locations of mask5 fit in
                    %make bigmask5 that spans whole cortex
                    
                    bigmask5 = zeros(l_bigy,l_bigx);
                    bigmask4 = bigmask5;
                    bigedge_mask = bigmask5;
                    
                    yind_lt_corner_mask5 = round(l_bigy*(maxcortY-maxrot_Y)/(maxcortY-mincortY));
                    xind_lt_corner_mask5 = round(l_bigx*(minrot_X-mincortX)/(maxcortX-mincortX));
                    
                    bigmask5(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(mask5);
                    bigmask4(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(mask4);
                    bigedge_mask(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(edge_mask);
                    
%                     yind_lt_corner_mask5 = round(l_bigy*(minrot_Y-mincortY)/(maxcortY-mincortY));
%                     xind_lt_corner_mask5 = round(l_bigx*(minrot_X-mincortX)/(maxcortX-mincortX));
%                     bigmask5(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(mask5);
%                     bigmask4(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(mask4);
%                     bigedge_mask(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(edge_mask);
                    
                    hold on;
                    h2 = imagesc(bigmask5,[0 3]);
                    colormap vga
                    set(h2,'AlphaData',.5);
                    
                    but = 1;
                    while but == 1 | but == 3
                        [xi,yi,but] = ginput(1);
                        if but==1 | but == 3
                            if round(yi) < (l_bigy-2) && round(yi) > 2 && round(xi) < (l_bigx-2) && round(xi) > 2
                                if but == 1
                                    bigmask4((round(yi)-2):(round(yi)+2),(round(xi)-2):(round(xi)+2)) = 1;
                                    %also change mask4
                                elseif but == 3
                                    bigmask4((round(yi)-2):(round(yi)+2),(round(xi)-2):(round(xi)+2)) = 0;
                                    %also change mask4
                                end
                                %imagesc(mask4+2*edge_mask);
                                set(h2,'CData',bigmask4+2*bigedge_mask)
                            end
                        end
                    end
                    allmasks{i,areanum,methodInd} = flipud(bigmask4);
                    %close figure
                    close(100)
                else
%                     %DISPLAY SAVED MASK
%                     figure(100)
%                     %imagesc(allmasks{i,areanum,methodInd})
%                     bigmask4 = flipud(allmasks{i,areanum,methodInd});
%                     h2 = imagesc(bigmask4,[0 3]);
%                     colormap vga
%                     set(h2,'AlphaData',.5);
%                     
%                     
%                     but = 1;
%                     while but == 1 | but == 3
%                         [xi,yi,but] = ginput(1);
%                         if but==1 | but == 3
%                             if round(yi) < (l_bigy-2) && round(yi) > 2 && round(xi) < (l_bigx-2) && round(xi) > 2
%                                 if but == 1
%                                     bigmask4((round(yi)-2):(round(yi)+2),(round(xi)-2):(round(xi)+2)) = 1;
%                                     %also change mask4
%                                 elseif but == 3
%                                     bigmask4((round(yi)-2):(round(yi)+2),(round(xi)-2):(round(xi)+2)) = 0;
%                                     %also change mask4
%                                 end
%                                 %imagesc(mask4+2*edge_mask);
%                                 set(h2,'CData',bigmask4)
%                             end
%                         end
%                     end
%                     allmasks{i,areanum,methodInd} = flipud(bigmask4);
%                     
%                     close(100)
                end
            
                fprintf(1, 'done\n')       
                
                bigmask4 = allmasks{i,areanum,methodInd};
                mask_perim = bwperim(bigmask4,8);
                mask_perim_nums = union(mask_perim_nums,[areanum]);
                mask_perims{areanum} = mask_perim;
                
                %-----------------------------------
                %12c. Plot whole cortex only when areanum=7
                if areanum==8
                    
                    %12ca. Plot cortex at different timepoints                    
                    figure(20+4*areanum)
                    
                    index=0;    
                    for timenow = timestart:timeinc:timestop
                        index=index+1;
                        
                        %16 subplots
                        subplot(4,4,index)
                        hold on
                        
                        %SURFACE PLOT
                        C=targetC(:,timenow);
                        CI = griddata(rot_targetL(:,1),rot_targetL(:,2),C,rot_XI,rot_YI);
                        
                        
                        surface('XData',rot_XI,'YData',rot_YI,'ZData',rot_ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6)
                        
                        %plot CENTRAL SULCUS FOR REFERENCE
                        plot3(rot_censul(:,1),rot_censul(:,2),rot_censul(:,3),'LineWidth',7,'color','g','LineStyle','.');
                        plot3(rot_censul(:,1),rot_censul(:,2),rot_censul(:,3),'Markersize',4,'color','g','Marker','.');
%                         plot3(diam1(:,1),diam1(:,2),diam1(:,3),'LineWidth',7,'color','r','LineStyle','.');
%                         plot3(diam1(:,1),diam1(:,2),diam1(:,3),'Markersize',4,'color','r','Marker','.');
                        
                        view(-90,90)
                        TR = timeRange_total{i};
                        title(['time = ',num2str(TR(1)+(TR(2)-TR(1))*(timenow-1)/(timePoints_total(i))),'ms'])
                        axis off
                        colorbar
                    end
                    
%                     %12d. Plot one whole averaged cortex over time period chosen
%                     figure(21+4*areanum)
%                     timenow=plotstart:plotstop;
%                     C = mean(targetC(:,timenow),2);
%                     hold on
%                     CI = griddata(targetL(:,1),targetL(:,2),C,XI,YI);
%                     surface('XData',XI,'YData',YI,'ZData',ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6)
%                     
%                     plot3(censul1(:,1),censul1(:,2),censul1(:,3),'LineWidth',7,'color','g','LineStyle','.');
%                     plot3(censul1(:,1),censul1(:,2),censul1(:,3),'Markersize',4,'color','g','Marker','.');
%                     plot3(censul2(:,1),censul2(:,2),censul2(:,3),'LineWidth',7,'color','g','LineStyle','.');
%                     plot3(censul2(:,1),censul2(:,2),censul2(:,3),'Markersize',4,'color','g','Marker','.');
%                     %plot3(diam1(:,1),diam1(:,2),diam1(:,3),'Markersize',4,'color','r','Marker','.');
%                     
%                     view(90,90)
%                     axis off
%                     colorbar
                end
                
                %-------------------------------------
                %12e. Make new edge and find new targetC and targetL
                %locations
                
                fprintf(1, 'create new mask... ')
                
                %CONVERT MASKS to LOCATION MASKS
                [mask_m,mask_n] = find(bigmask4>0);
                for mask_ind = 1:length(mask_m)
                    rot_new_edge(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
                end
                
                [mask_m,mask_n] = find(mask_perim>0);
                for mask_ind = 1:length(mask_m)
                    rot_new_perim(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
                end
                
                regions{i,areanum,methodInd} = rot_new_perim;
                alledges{i,areanum,methodInd} = rot_new_edge;
                fprintf(1, 'done\n')
                
                
                %GET NEW STRENGTHS AND BASELINE STRENGTHS
                fprintf(1, 'get new strengths... ')
                %[newtargetL,newtargetC]=find_new_cortex_strengths_AC(top_rot_cortexL,top_cortexC,rot_new_edge,timePoints,0);
                [newtargetL,newtargetC]=find_toplayer_cortex_strengths_AC(top_rot_cortexL,top_cortexC,rot_new_edge,rot_edgeL,timePoints,0);
                fprintf(1, 'done\n')
                
                fprintf(1, 'get new baseline strengths...')
                [newtargetL_BL,newtargetC_BL]=find_new_cortex_strengths_AC(top_rot_cortexL,top_cortexC_BL,rot_new_edge,timePoints_BL,0);
                fprintf(1, 'done\n')
                
                
                
                
                [m,n] = size(newtargetC_BL)
                newlongtargetC_BL = reshape(newtargetC_BL,m*n,1);
                newsorttargetC_BL = sort(newlongtargetC_BL);
                newstdtargetC_BL = std(newsorttargetC_BL);
                
                [m,n] = size(newtargetC)
                newlongtargetC = reshape(newtargetC,m*n,1);
                newsorttargetC = sort(newlongtargetC);
                
                maxtargetC = max(max(newtargetC(:,1:plotstop)));
                maxvalues(areanum) = maxtargetC;
                                
                newThreshold = max(newsorttargetC_BL)
                newThreshold = .11*maxtargetC
                newThresholds(areanum) = newThreshold;
                
                thresholds = newThreshold;
                
                PercentThresholdofMax = newThreshold/maxtargetC
                percentmaxthresholds(areanum) = PercentThresholdofMax;
                
%                 figure(90+i)
%                 plot(newsorttargetC_BL,'Color','g','Linewidth',2)
%                 hold on;
%                 plot(newsorttargetC,'Color','b','Linewidth',2)
%                 line([1 length(newsorttargetC)],[newThreshold newThreshold],'Color','r','Linewidth',2)
%                             
                
                
                %12f. Plot one whole averaged cortex over time period chosen
                figure(22+4*areanum+1)
                timenow=1:plotstop;
                C = mean(newtargetC(:,timenow),2);
                maxmeanC = max(max(C));
                fprintf('max mean strength = %d\n',maxmeanC)
                hold on
                
                
                %unrotate points back, so that can create a .pom file with
                %new locations
                unRx = [1 0 0; 0 cos(-alpha) sin(-alpha); 0 -sin(-alpha) cos(-alpha)];
                unRy = [cos(-beta) 0 -sin(-beta); 0 1 0; sin(-beta) 0 cos(-beta)];
                unRz = [cos(-gamma) sin(-gamma) 0; -sin(-gamma) cos(-gamma) 0; 0 0 1];
                unRx = unRx*unRy*unRz;
                unrot_newtargetL = (unRx*newtargetL')';
                
                %create new .pom file of new locations
                write_new_pom_file(areanames,areanum,Base_dir,subjectName,unrot_newtargetL);
                
                
                %find limits of new target locations
                newminrot_X=min(newtargetL(:,1));
                newmaxrot_X=max(newtargetL(:,1));
                newminrot_Y=min(newtargetL(:,2));
                newmaxrot_Y=max(newtargetL(:,2));
                newminrot_Z=min(newtargetL(:,3));
                newmaxrot_Z=max(newtargetL(:,3));
                newrot_x=[newminrot_X:samp_inc:newmaxrot_X];
                newrot_y=[newminrot_Y:samp_inc:newmaxrot_Y];
                [newrot_XI,newrot_YI] = meshgrid(newrot_x,newrot_y); 
                
                CI = griddata(newtargetL(:,1),newtargetL(:,2),C,newrot_XI,newrot_YI);
                newrot_ZI = griddata(newtargetL(:,1),newtargetL(:,2),newtargetL(:,3),newrot_XI,newrot_YI);
                
                %SURFACE PLOT OVERLAYED WITH MASK PERIM
                surface('XData',newrot_XI,'YData',newrot_YI,'ZData',newrot_ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6)
                
                plot3(rot_new_perim(:,1),rot_new_perim(:,2),rot_new_perim(:,3),'LineWidth',7,'color','r','LineStyle','.')
                
                %scatter3(rot_censul(:,1),rot_censul(:,2),rot_censul(:,3),3,'g');
                                
                view(0,90)
                axis off
                axis equal
                colorbar
                v_caxis = caxis;
                
                %-------------------------
                targetCI = [];
%                 for timenow = 1:timePoints
%                     C=newtargetC(:,timenow);
%                     CI = griddata(newtargetL(:,1),newtargetL(:,2),C,XI,YI);
%                                         
%                     [m,n] = size(CI);
%                     dummyCI = reshape(CI,m*n,1);
%                     targetCI(:,timenow) = dummyCI(~isnan(dummyCI));
%                 end
                
                targetCI = newtargetC;

                %---------------------------------------------------
                %13a. Figure out statistical threshold of activation
                
%                 [m,n] = size(targetCI);
%                 longtargetCI = reshape(targetCI,m*n,1);
%                 sorttargetCI = sort(longtargetCI);
%                 
%                 for sigtime = 1:200:length(sorttargetCI)
%                     T1 = length(sorttargetCI)-sigtime;
%                     T2 = sigtime;
%                     group1=repmat(V1,T1,1);
%                     group2=repmat(V2,T2,1);
%                     
%                     cellGroup=(cellstr([group1;group2]))';
%                     
%                     %use ANOVA to calculate significance
%                     p_thr(length(sorttargetCI)-sigtime+1) = anova1(sorttargetCI,cellGroup,'off');
%                     
%                 end
%                 %f1 = gcf;
%                 
%                 figure
%                 plot(sorttargetCI)
%                 totaltime = timeRange(2)-timeRange(1);
%                 times1 = timeRange(1):totaltime/(timePoints-1):timeRange(1)+50;
%                 
%                 [pmin,pind] = min(p_thr);
%                 if ~isempty(pind) && pind<length(sorttargetCI)
%                     thr_act = sorttargetCI(pind)
%                     line([pind pind],[0 max(sorttargetCI)],'Color','r','Linewidth',2)
%                     %figure(f1)
%                 end
                
                               
                %--------------------------------------------------------
                %13. Statistical determination of significant differences
                %from baseline over range of time window
                
%                 %number of points
%                 %[locationNo_BL,TimeNo_BL]=size(targetC_BL);
%                 [locationNo,TimeNo]=size(targetC);
%                 
%                 %loop through all points and determine whether significant
%                 % difference from baseline over the range of time window
%                 siglocationInd = [];
%                 for locationInd=1:locationNo
%                     X=[targetC_BL(locationInd,:),targetC(locationInd,:)];
%                     group1=repmat(V1,TimeNo_BL,1);
%                     group2=repmat(V2,TimeNo,1);
%                     
%                     cellGroup=(cellstr([group1;group2]))';
%                     
%                     %use ANOVA to calculate significance
%                     p(methodInd,locationInd) = anova1(X,cellGroup,'off');
%                     
%                     %if significance below value, add to count
%                     if p(methodInd,locationInd)<0.05
%                         ActCount=ActCount+1;
%                         siglocationInd = [siglocationInd,locationInd];
%                     end    
%                 end %end of locationInd
%                 
%                 ActiveRatio(k,i)=ActCount/locationNo;
%                 %end loop through points to determine significance
                %----------------------------------------------------------
                
                %14. Calculate baseline average
                
                
                
                %----------------------------------------------------------
                %15. Calculate threshold areas inside each region above baseline average and
                %plot area vs time and avg strengths vs time
                %look for maximum value in time window up to 0ms
                %maxvalue{i,areanum,methodInd} = max(max(targetCI(:,1:plotstop)));
                                
                %thresholds are a percentage of maximum value inside that
                %area   %DOES NOT USE FIRST SET OF DATA FOR MAX VALUE TO STANDARDIZE
                %ACROSS DIFFERENT SETS OF DATA FOR THE SAME PERSON

                %thresholds = [.85; .80; .75; .70; .65; .60; .50; .40; .20]*maxvalue{i,areanum,methodInd};
                                
                totaltime = timeRange(2)-timeRange(1);
                times1 = timeRange(1):totaltime/(timePoints-1):timeRange(1)+50;
                
                fprintf(1, 'calculating indices... ')
                for oind = 1:length(thresholds)
                    %number of locations in that area above threshold value
                    overthresholds{oind,i,areanum,methodInd} = sum(targetCI>thresholds(oind),1);
                    
                    %percentage of active area that is over threshold
                    %totalarea = size(targetCI,1); %use imageprocessing way to do this?
                    totalarea = sum(targetCI>0,1);
                    totalarea = totalarea(1);
                    percentactiveareas{oind,i,areanum,methodInd} = overthresholds{oind,i,areanum,methodInd}/totalarea;
                                        
                    %total strengths of active areas over threshold
                    [t_ind1,t_ind2,v] = find(targetCI>thresholds(oind));
                    totalsum = zeros(1,timePoints);
                    for ind = 1:length(t_ind1)
                        totalsum(1,t_ind2(ind)) = totalsum(t_ind2(ind))+targetCI(t_ind1(ind),t_ind2(ind));
                    end
                    totalstrengths{oind,i,areanum,methodInd} = totalsum;
                    %totalstrengths{oind,i,areanum,methodInd} = sum(targetCI(targetCI>thresholds(oind)),1);
                    
                    %average strength of active areas over threshold
                    %dummy_avgstrengths = totalstrengths{oind,i,areanum,methodInd}./overthresholds{oind,i,areanum,methodInd}./maxvalue{i,areanum,methodInd};
                    dummy_avgstrengths = totalstrengths{oind,i,areanum,methodInd}./overthresholds{oind,i,areanum,methodInd};
                    dummy_avgstrengths(isnan(dummy_avgstrengths)) = 0;
                    avgstrengths{oind,i,areanum,methodInd} = dummy_avgstrengths;
                    
                    %hemisphere laterality index is weighted measure of strength and location
                    HLIs{oind,i,areanum,methodInd} = totalstrengths{oind,i,areanum,methodInd}.*overthresholds{oind,i,areanum,methodInd};

                    %"center of activity" is strength weighted locations, averaged together 
                    %strengths = overthresholds{oind,i,areanum,methodInd}
                    %locations = targetL;
                    [m_target,n_target] = size(targetCI);
                    
                    for col = 1:n_target
                        S_col = targetCI(:,col);
                        sum_S_col = sum(S_col,1);
                        S_col_rep = repmat(S_col,1,3);
                        wt_loc = newtargetL.*S_col_rep;
                        mean_loc(:,col) = (1/sum_S_col)*sum(wt_loc,1)';
                    end
                    single_mean_loc = mean(mean_loc,2);
                    Subject.single_mean_loc = single_mean_loc;
                    
                    %plot where center of activity occurs
                    figure(101)
                    hold on;
                    plot3(rot_new_perim(:,1),rot_new_perim(:,2),rot_new_perim(:,3),'LineWidth',7,'color','r','LineStyle','.')
                    %scatter3(mean_loc(1,:),mean_loc(2,:),mean_loc(3,:),10,0:2/(length(mean_loc(1,:))-1):2);
                    plot3(single_mean_loc(1),single_mean_loc(2),single_mean_loc(3),'Markersize',30,'color','m','LineStyle','.')
                    %scatter3(rot_targetL(:,1),rot_targetL(:,2),rot_targetL(:,3),7,'m');
                    %scatter3(newtargetL(:,1),newtargetL(:,2),newtargetL(:,3),7,'y');
                    midline3_x = P(1)*newrot_y+P(2);
                    plot(midline3_x,newrot_y,'color','r','Linewidth',2,'Linestyle','--')
                    
                    %find distance to midline
                    u_vector = [single_mean_loc(1)-midline3_x(1),single_mean_loc(2)-newrot_y(1)];
                    v_vector = [midline3_x(end)-midline3_x(1),newrot_y(end)-newrot_y(1)];
                    
                    proj_v_u = dot(v_vector,u_vector)/(dot(v_vector,v_vector))*v_vector;
                    
                    orth_midline_pt = [midline3_x(1) newrot_y(1)] + proj_v_u;

                    line([single_mean_loc(1) orth_midline_pt(1)],[single_mean_loc(2) orth_midline_pt(2)],'Color','b')
                    dist_to_midline = sqrt((single_mean_loc(1)-orth_midline_pt(1))^2+(single_mean_loc(2)-orth_midline_pt(2))^2)
                    
                    view(0,90)
                    axis on
                    axis equal
                    %title(['distance to midline = ',num2str(dist_to_midline)])
                                        
                    
                    %overlay region of interest onto picture of cortex
%                     imA = imread('WScortex.jpg');
%                     imB = im2bw(imA,.95);
%                     se = strel('square',5);
%                     imC = imclose(imB,se);
%                     imD = ~imC;
%                     imE = imfill(imD,'holes');
%                     [im_i,im_j] = find(imE>0);
%                     imF = imE(min(im_i):max(im_i),min(im_j):max(im_j));
%                     imG = bwperim(imF,8);
%                     
%                     imA_resize = imA(min(im_i):max(im_i),min(im_j):max(im_j));
%                     %fit cortex picture to cortex pt locations
%                     
%                     rot_cortexL = (Rx*cortexL')';
%                     
%                     mincortX=min(rot_cortexL(:,1));
%                     maxcortX=max(rot_cortexL(:,1));
%                     mincortY=min(rot_cortexL(:,2));
%                     maxcortY=max(rot_cortexL(:,2));
%                     
%                     [imG_m,imG_n] = size(imG);
%                     imA_resize2 = imresize(imA_resize,[imG_m,imG_n],'bilinear');
%                     
%                     samp_incX = (maxcortX-mincortX)/(imG_n-1);
%                     samp_incY = (maxcortY-mincortY)/(imG_m-1);
%                     
%                     x_locs = mincortX:samp_incX:maxcortX;
%                     y_locs = maxcortY:-samp_incY:mincortY;
%                     
%                     [imG_pts_i,imG_pts_j] = find(imG>0);
%                     for ind_cort=1:5:length(imG_pts_i)
%                         loc_cortex(ind_cort,:) = [x_locs(imG_pts_j(ind_cort)),y_locs(imG_pts_i(ind_cort))];
%                     end
%                     loc_cortex = loc_cortex(1:5:end,:);                
%                     
%           
%                     %rotate everything else too
%                     rot_new_perim = (Rx*new_perim')';
%                     rot_censul = (Rx*censul')';
%                                          
                    figure(102)
                    imshow([mincortX maxcortX],[maxcortY mincortY],imA_resize2)
                    axis on
                    axis([mincortX maxcortX mincortY maxcortY])
                    hold on;
                    %h3 = surface('XData',newrot_XI,'YData',newrot_YI,'ZData',newrot_ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6)
                    %set(h3,'AlphaData',.6);
                    
                    %scatter3(mean_loc(1,:),mean_loc(2,:),mean_loc(3,:),10,'g');
                    %scatter3(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),7,'b');
                    %scatter3(loc_cortex(:,1),loc_cortex(:,2),zeros(length(loc_cortex(:,1)),1),7,'m')
                    plot(rot_new_perim(:,1),rot_new_perim(:,2),'LineWidth',7,'color','r','LineStyle','.')
                    
%                     for i_mask = 1:length(mask_perim_nums)
%                         dum_region = regions{i,i_mask,methodInd};
%                         plot(dum_region(:,1),dum_region(:,2),'LineWidth',7,'color','r','LineStyle','.')
%                     end
                    
                    %scatter(rot_censul(:,1),rot_censul(:,2),15,'g');
                    
                    view(0,90)
                    axis equal
                    axis xy
                    caxis(v_caxis);
%                     colorbar
                    
                    
                    
                    %set NaN entries to 0 if no values above threshold
                    dummy_avg = avgstrengths{oind,i,areanum,methodInd};
                    dummy_avg(isnan(dummy_avg)) = 0;
                    
                    %subtract baseline from strengths
                    %basestrengths{oind,i,areanum,methodInd} = mean(dummy_avg(1:length(times1)));
                    %avgstrengths{oind,i,areanum,methodInd} = (dummy_avg - basestrengths{oind,i,areanum,methodInd});

                    %maxv = maximum value seen
                    maxv = maxtargetC;
                    SIs{oind,i,areanum,methodInd} = max(avgstrengths{oind,i,areanum,methodInd})/maxv;

                    %active area index
                    AAIs{oind,i,areanum,methodInd} = max(percentactiveareas{oind,i,areanum,methodInd});
                    
                    %subtract baseline from areas
                    %dummy_area = percentactiveareas{oind,i,areanum,methodInd};
                    %baseareas{oind,i,areanum,methodInd} = mean(dummy_area(1:length(times1)));
                    %percentactiveareas{oind,i,areanum,methodInd} = (dummy_area - baseareas{oind,i,areanum,methodInd});

                end
                fprintf(1, 'done\n')
                
                
                if areanum < 7
                    %regions graph of active areas as percentage of total area vs time
                    if areanum<=4
                        figure(i+6)
                        subplot(3,1,areanum)    %SMA,PM,M1,S1
                    elseif areanum>4 && areanum <=6
                        figure(i+9)             %start at figure 10, end at 12
                        subplot(2,1,areanum-4)  %right and left hemisphere
                    end
                    hold on;
                    set(get(gca,'YLabel'),'Rotation',0.0);
                    
                    fprintf(1, 'graphing active area... ')
                    for oind = 1:length(thresholds)
                        
                        if oind == 1
                            color = [255 255 102]/255;
                            %color = [125 190 102]/255;
                        elseif oind == 2
                            color = [243 249 102]/255;
                        elseif oind == 3
                            color = [226 241 102]/255;
                        elseif oind == 4
                            color = [206 231 102]/255;
                        elseif oind == 5
                            color = [170 212 102]/255;
                        elseif oind == 6
                            color = [125 190 102]/255;
                        elseif oind == 7
                            color = [97 176 102]/255;
                        elseif oind == 8
                            color = [40 148 102]/255;
                        elseif oind == 9
                            color = [0 128 102]/255;
                        end
                        
                        area(timeRange(1):totaltime/(timePoints-1):timeRange(2),percentactiveareas{length(thresholds)-oind+1,i,areanum,methodInd},'facecolor',color);
                        times1 = timeRange(1):totaltime/(timePoints-1):0;
                        n_dum = percentactiveareas{length(thresholds)-oind+1,i,areanum,methodInd};
                        meanareatimes{oind,i,areanum,methodInd} = dot(times1,n_dum(1:length(times1)))/sum(n_dum(1:length(times1)));
                        sdareatimes{oind,i,areanum,methodInd} = sqrt((sum(n_dum(1:length(times1)).*(times1-meanareatimes{oind,i,areanum,methodInd}).^2))/(sum(n_dum(1:length(times1)))-1));
                        
                        %line([(meanareatimes{oind,i,areanum,methodInd}-sdareatimes{oind,i,areanum,methodInd}) (meanareatimes{oind,i,areanum,methodInd}+sdareatimes{oind,i,areanum,methodInd})],[max(n_dum) max(n_dum)],'Color',color,'linewidth',2)
                        
                        %figure out where significant activation begins
                        X=n_dum;
                        %T1 = times for nonactivation
                        %T2 = times for activation
                        %                     fprintf(1, 'calculating significant activation time onset(%d)... ',oind)
                        %                     for sigtime = 1:round(length(X)/200):length(X)
                        %                         T1 = length(X)-sigtime;
                        %                         T2 = sigtime;
                        %                         group1=repmat(V1,T1,1);
                        %                         group2=repmat(V2,T2,1);
                        %                     
                        %                         cellGroup=(cellstr([group1;group2]))';
                        %                         
                        %                         %use ANOVA to calculate significance
                        %                         p(length(X)-sigtime+1) = anova1(X,cellGroup,'off');
                        %                         
                        %                     end
                        %                     fprintf(1,'done\n')
                        %                     
                        %                     %f1 = gcf;
                        %                     %figure(2)
                        %                     %plot(p')
                        %                     [pmin,pind] = min(p);
                        %                     if ~isempty(pind) && pind<length(times1)
                        %                         areaacttimes{length(thresholds)-oind+1,i,areanum,methodInd} = times1(pind);
                        %                         line([times1(pind) times1(pind)],[0 max(n_dum)],'Color','r','Linewidth',2)
                        %                         %figure(f1)
                        %                     end
                    end
                    fprintf(1, 'done\n')
                    
                    %generate nice graph for brad
                    %                 pa5 = percentactivearea5{i,areanum,methodInd};
                    %                 if areanum==1
                    %                     pa5(1:50) = 0;
                    %                 elseif areanum==2
                    %                 elseif areanum==3
                    %                     pa5(1:65) = 0;
                    %                 elseif areanum==4
                    %                     pa5(1:77) = 0;
                    %                 end 
                    %                 a5 = area(timeRange(1):totaltime/(timePoints-1):timeRange(2),pa5,'facecolor',[170 212 102]/255);    
                    
                    if areanum==1
                        ylabel('AAI_S_M_A          ','Fontsize',14)
                    elseif areanum==2
                        ylabel('AAI_P_M            ','Fontsize',14)   
                    elseif areanum==3
                        ylabel('AAI_M_1            ','Fontsize',14)
                    elseif areanum==4
                        ylabel('AAI_S_1            ','Fontsize',14)
                        xlabel('time relative to movement onset (ms)','Fontsize',14)
                    elseif areanum==5
                        ylabel('AAI_R_H            ','Fontsize',14)
                    elseif areanum==6
                        ylabel('AAI_L_H            ','Fontsize',14)
                        xlabel('time relative to movement onset (ms)','Fontsize',14);    
                    end
                    
                    %title('active area as percentage of region of interest above threshold')
                    %legend('80%','60%','40%','20%',-1)
                    
                    %regions graph of avg strengths vs time
                    
                    if areanum<=4
                        figure(i)
                        subplot(4,1,areanum)    %SMA,PM,M1,S1
                    else
                        figure(i+3)             %start at figure 4, end at 6
                        subplot(2,1,areanum-4)  %right and left hemisphere
                    end
                    hold on;
                    set(get(gca,'YLabel'),'Rotation',0.0);
                    
                    fprintf(1, 'graphing strengths... ')
                    for oind = 1:length(thresholds)
                        n_avgstrengths{oind,i,areanum,methodInd} = avgstrengths{oind,i,areanum,methodInd}/maxv;
                        
                        if oind == 1
                            color = [0 128 102]/255;
                        elseif oind == 2
                            color = [40 148 102]/255;
                        elseif oind == 3
                            color = [97 176 102]/255;
                        elseif oind == 4
                            color = [125 190 102]/255;
                        elseif oind == 5
                            color = [170 212 102]/255;
                        elseif oind == 6
                            color = [206 231 102]/255;
                        elseif oind == 7
                            color = [226 241 102]/255;
                        elseif oind == 8
                            color = [243 249 102]/255;
                        elseif oind == 9
                            color = [255 255 102]/255;
                        end
                        
                        area(timeRange(1):totaltime/(timePoints-1):timeRange(2),n_avgstrengths{oind,i,areanum,methodInd},'facecolor',color);
                        
                        %graph temporally where mean activity occurs at what timepoint/timerange
                        
                        times1 = timeRange(1):totaltime/(timePoints-1):0;
                        n_dum = n_avgstrengths{oind,i,areanum,methodInd};
                        meantimes{oind,i,areanum,methodInd} = dot(times1,n_dum(1:length(times1)))/sum(n_dum(1:length(times1)));
                        sdtimes{oind,i,areanum,methodInd} = sqrt((sum(n_dum(1:length(times1)).*(times1-meantimes{oind,i,areanum,methodInd}).^2))/(sum(n_dum(1:length(times1)))-1));
                        
                        %line([(meantimes{oind,i,areanum,methodInd}-sdtimes{oind,i,areanum,methodInd}) (meantimes{oind,i,areanum,methodInd}+sdtimes{oind,i,areanum,methodInd})],[max(n_dum) max(n_dum)],'Color',color,'linewidth',2)
                        
                        %figure out where significant activation begins
                        X=n_dum;
                        %T1 = times for nonactivation
                        %T2 = times for activation
                        
                        %                     disp 'calculating significant activation time onset'
                        %                     for sigtime = 1:round(length(X)/200):length(X)
                        %                         T1 = length(X)-sigtime;
                        %                         T2 = sigtime;
                        %                         group1=repmat(V1,T1,1);
                        %                         group2=repmat(V2,T2,1);
                        %                     
                        %                         cellGroup=(cellstr([group1;group2]))';
                        %                         
                        %                         %use ANOVA to calculate significance
                        %                         p(length(X)-sigtime+1) = anova1(X,cellGroup,'off');
                        %                         
                        %                     end
                        %                     disp 'done'
                        %                     %f1 = gcf;
                        %                     %figure(2)
                        %                     %plot(p')
                        %                     [pmin,pind] = min(p);
                        %                     if ~isempty(pind) && pind<length(times1)
                        %                         strengthacttimes{oind,i,areanum,methodInd} = times1(pind);
                        %                         line([times1(pind) times1(pind)],[0 max(n_dum)],'Color','r','Linewidth',2)
                        %                         %figure(f1)
                        %                     end
                        
                    end    
                    fprintf(1, 'done\n')
                    
                    %title('avg strength of active area above threshold')
                    %legend([a1,a2,a3,a4,a5,a6,a7,a8,a9],'85%','80%','75%','70%','65%','60%','50%','40%','20%',-1)
                    
                    if areanum==1
                        ylabel('SI_S_M_A          ','Fontsize',14)
                    elseif areanum==2
                        ylabel('SI_P_M            ','Fontsize',14)
                    elseif areanum==3
                        ylabel('SI_M_1            ','Fontsize',14)
                    elseif areanum==4
                        ylabel('SI_S_1            ','Fontsize',14)
                        xlabel('time relative to movement onset (ms)','Fontsize',14);
                    elseif areanum==5
                        ylabel('SI_R_H            ','Fontsize',14)
                    elseif areanum==6
                        ylabel('SI_L_H            ','Fontsize',14)
                        xlabel('time relative to movement onset (ms)','Fontsize',14);
                    end
                    
                    %end calculate thresholded areas
                    %----------------------------------------------------------
                    
                    
                    %----------------------------------------------------------
                    %graph average strength of significant areas inside region versus time
                    %                 %avgstrength{areanum,methodInd} = mean(targetC(siglocationInd',:),1)
                    %                 %avgstrength_BL{areanum,methodInd} = mean(targetC_BL(siglocationInd',:),1)
                    %                 avgstrength{areanum,methodInd} = mean(targetC,1)
                    %                 avgstrength_BL{areanum,methodInd} = mean(targetC_BL,1)
                    %                 
                    %                 %regions graph
                    %                 figure(3)
                    %                 hold on;
                    %                 subplot(4,1,areanum)
                    %                 plot(-300:350/(timePoints-1):50,avgstrength{areanum,methodInd},'r')
                    %                 
                    %                 %baseline graph
                    %                 figure(6)
                    %                 hold on;
                    %                 subplot(4,1,areanum)
                    %                 plot(-1900:200/(timePoints_BL-1):-1700,avgstrength_BL{areanum,methodInd},'b')
                    %end graph avg strength
                    %----------------------------------------------------------
                    
                end  
            end %end of methodInd
            %--------------------------------------------------------------
        
        end %end of areanum
        %------------------------------------------------------------------
            
    end % end of task
    %----------------------------------------------------------------------
    
end % end of subject Number

%display active ratio
%ActiveRatio

beep

%random things to display

%%calculate Hemisphere Location Index
if areanum < 7
    for i=1:1
        size_HLIs = size(HLIs);
        if length(size_HLIs)>2
            if size_HLIs(3)>=6
                figure(i+12)
                hold on
                for oind = 1:length(thresholds)
                    
                    HLIcontra = HLIs{oind,i,5,1};   %right(5) is contra 
                    HLIipsi = HLIs{oind,i,6,1};     %left(6) is ipsi
                    HLIindex1 = (HLIcontra-HLIipsi)./(HLIcontra+HLIipsi);
                    HLIindex1(isnan(HLIindex1))=0;
                    
                    HLIindices{oind,i,methodInd} = HLIindex1;
                    meanHLIindices{oind,i,methodInd} = mean(HLIindex1(1:plotstop));
                    if oind == 1
                        color = [0 128 102]/255;
                    elseif oind == 2
                        color = [40 148 102]/255;
                    elseif oind == 3
                        color = [97 176 102]/255;
                    elseif oind == 4
                        color = [125 190 102]/255;
                    elseif oind == 5
                        color = [170 212 102]/255;
                    elseif oind == 6
                        color = [206 231 102]/255;
                    elseif oind == 7
                        color = [226 241 102]/255;
                    elseif oind == 8
                        color = [243 249 102]/255;
                    elseif oind == 9
                        color = [255 255 102]/255;
                    end
                    
                    a1 = area(timeRange(1):totaltime/(timePoints-1):timeRange(2),HLIindex1,'facecolor',color);
                end
            end
        end
    end
end

Subject.maxvalues = maxvalues;
Subject.newThresholds = newThresholds;
Subject.percentmaxthresholds = percentmaxthresholds;
Subject.SIs = SIs;
Subject.AAIs = AAIs;

if areanum < 7
    size_HLIs = size(HLIs);
    if length(size_HLIs)>2
        if size_HLIs(3)>=6
            Subject.HLIindices = HLIindices;
            Subject.meanHLIindices = meanHLIindices;
        end
    end
end
Subject.allmasks = allmasks;
Subject.mask_perims = mask_perims;
Subject.mask_perim_nums = mask_perim_nums;
Subject.regions = regions;










% 

% for i=1:1   %Strengths graph
%     figure(i)
%     for j=1:4
%         subplot(4,1,j)
%         
%         timeRange = timeRange_total{mod(i-1,3)+1};
%         axis([timeRange(1) timeRange(2) 0 2.5]);
%         %axis([-100 100 0 .01]);
%         grid on;
%     end
% end
% for i=7:7   %Areas graph
%     figure(i)
%     for j=1:4
%         subplot(4,1,j)
%         
%         timeRange = timeRange_total{mod(i-1,3)+1};
%         axis([timeRange(1) timeRange(2) 0 .025]);
%         %axis([-100 100 0 .01]);
%         grid on;
%     end
% end

% figure  %plot all regions on graph
% hold on;
% for j=1:4
%     region = regions{1,j,1};
%     edge = alledges{1,j,1};
%     plot3(edge(:,1),edge(:,2),edge(:,3),'LineWidth',7,'color','b','LineStyle','.')
%     plot3(region(:,1),region(:,2),region(:,3),'LineWidth',7,'color','r','LineStyle','.')
%     
% end
% view(90,90)




%SUBFUNCTIONS--------------------------------------------------------------

function [contrib,targetL,targetC] = find_cortex_strengths_AC(cdr_file_name,cortexL,edge,TimePoints,plotCort)

contrib=[];
fpos=0;

for i=1:TimePoints
    if i==1
        [contrib_tmp,Ccount,CNR]=read_Curry_file4_AC([cdr_file_name(1:end-4),'_strength.cdr'],'STRENGTH',1,0);
        % disp 'done1'
    elseif i==TimePoints
        [contrib_tmp,Ccount,CNR]=read_Curry_file4_AC([cdr_file_name(1:end-4),'_strength.cdr'],'STRENGTH',1,1);
        % disp 'done2'
    else
        [contrib_tmp,Ccount,CNR]=read_Curry_file4_AC([cdr_file_name(1:end-4),'_strength.cdr'],'STRENGTH',0,1);
        % disp 'done3'
    end
    contrib=[contrib,contrib_tmp];
end

%--------Select the points inside the edge---------------
Tind=[];
[n,tmp]=size(edge);
[m,tmp]=size(cortexL);

for i=1:n
    dis= sqrt(sum((cortexL-repmat(edge(i,:),size(cortexL,1),1)).^2,2));
    
    ind = find(dis<2);      %if want to look around pts
    %[yind,ind] = min(dis); %if want to pick exact pt that corresponds to selected pt in Curry
    
    if dis(ind)<5
        if i==1
            Tind = ind;
        end
        Tind = [Tind;ind];
    end
end
Tind = unique(Tind);

targetL=cortexL(Tind,:);
targetC=contrib(Tind,:);

if plotCort

    %plot3(cortexL(:,1),cortexL(:,2),cortexL(:,3),'b.');
	%hold on
	%plot3(targetL(:,1),targetL(:,2),targetL(:,3),'m+');
    
    figure
    tpind = 0;
    for tp=1:10:TimePoints
        tpind = tpind+1;
        subplot(4,ceil(TimePoints/(4*10)),tpind)
        %%%%plot3(cortexL(:,1),cortexL(:,2),cortexL(:,3),'b.');
        scatter3(targetL(:,1),targetL(:,2),targetL(:,3),20,targetC(:,tp),'filled');

        %-----------------------------------
        %set axis and colorbar properties for each person
        axis([-80 80 -75 75])  %for CM
        %axis([-100 100 -70 20])  %for BHrobot
        %caxis([0 2.0]);         %for BHrobot
        %------------------------------------
    
        axis equal
        hold on;
        
        colorbar
    end
    
    meantargetC = mean(targetC,2);
    
%     if TimePoints==52   %baseline graph, colormap stuff down below doesn't work?
%         map = colormap;
%         figure
%         colormap(map)
%     else
%         figure
%     end

    figure
    
    scatter3(targetL(:,1),targetL(:,2),targetL(:,3),20,meantargetC,'filled');
    axis([-80 80 -75 75])  %for CM
    %axis([-100 100 -70 20])  %for BHrobot
    axis equal
    
    %---------------------------
    %caxis([0 2.0]); %for BHrobot
    %---------------------------
    
    colorbar
    title('mean activity in time window')
end


function [targetL,targetC] = find_new_cortex_strengths_AC(cortexL,contrib,new_edge,TimePoints,plotCort)

fpos=0;

%--------Select the points inside the edge---------------
Tind=[];
[n,tmp]=size(new_edge);
[m,tmp]=size(cortexL);

for i=1:n
    dis= sqrt(sum((cortexL-repmat(new_edge(i,:),size(cortexL,1),1)).^2,2));
    
    %ind = find(dis<5);
    [yind,ind] = min(dis);
    if dis(ind)<5 
        if i==1
            Tind = ind;
        end
        Tind = [Tind;ind];
    end
end
Tind = unique(Tind);

targetL=cortexL(Tind,:);
targetC=contrib(Tind,:);

if plotCort

    %plot3(cortexL(:,1),cortexL(:,2),cortexL(:,3),'b.');
	%hold on
	%plot3(targetL(:,1),targetL(:,2),targetL(:,3),'m+');
    
    figure
    tpind = 0;
    for tp=1:10:TimePoints
        tpind = tpind+1;
        subplot(4,ceil(TimePoints/(4*10)),tpind)
        %%%%plot3(cortexL(:,1),cortexL(:,2),cortexL(:,3),'b.');
        scatter3(targetL(:,1),targetL(:,2),targetL(:,3),20,targetC(:,tp),'filled');

        %-----------------------------------
        %set axis and colorbar properties for each person
        axis([-80 80 -75 75])  %for CM
        %axis([-100 100 -70 20])  %for BHrobot
        %caxis([0 2.0]);         %for BHrobot
        %------------------------------------
    
        axis equal
        hold on;
        
        colorbar
    end
    
    meantargetC = mean(targetC,2);
    
%     if TimePoints==52   %baseline graph, colormap stuff down below doesn't work?
%         map = colormap;
%         figure
%         colormap(map)
%     else
%         figure
%     end

    figure
    
    scatter3(targetL(:,1),targetL(:,2),targetL(:,3),20,meantargetC,'filled');
    axis([-80 80 -75 75])  %for CM
    %axis([-100 100 -70 20])  %for BHrobot
    axis equal
    
    %---------------------------
    %caxis([0 2.0]); %for BHrobot
    %---------------------------
    
    colorbar
    title('mean activity in time window')
end

function [targetL,targetC] = find_toplayer_cortex_strengths_AC(cortexL,contrib,new_edge,old_edge,TimePoints,plotCort)

fpos=0;

%--------Select the points inside the edge---------------
Tind=[];
[n,tmp]=size(new_edge);
[m,tmp]=size(cortexL);

for i=1:n
    dis2total_edge =  sqrt(sum((old_edge-repmat(new_edge(i,:),size(old_edge,1),1)).^2,2));
    
    if min(dis2total_edge)<10
        dis2cortex= sqrt(sum((cortexL-repmat(new_edge(i,:),size(cortexL,1),1)).^2,2));
        
        [yind,ind] = min(dis2cortex);
        if dis2cortex(ind)<5
            if i==1
                Tind = ind;
            end
            Tind = [Tind;ind];
        end
    end
end
Tind = unique(Tind);

targetL=cortexL(Tind,:);
targetC=contrib(Tind,:);


function write_new_pom_file(areanames,areanum,Base_dir,subjectName,unrot_newtargetL)

num_locs = length(unrot_newtargetL(:,1));

pom_name = areanames{areanum};
fid1 = fopen([Base_dir,subjectName,'\Data\MRI\','new_',pom_name,'.pom'],'w');

fprintf(fid1,'%s\n','POINT_KEYWORDS START	# Do not edit!');
fprintf(fid1,'%s\n','POINT_KEY_LOCATIONS  = LOCATION_LIST');
fprintf(fid1,'%s\n','POINT_KEY_NORMALS    = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_CONTRIB    = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_FLAGS      = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_STRENGTHS  = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_ERRORS     = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_DEVIATIONS = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_FIELDS     = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_MGFP       = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_ADDITIVE   = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_MULTIPLICATIVE = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_PCA        = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_COLORIND   = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_CHARTRAFO  = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_NUMBERS    = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_NEIGHBORS  = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_TRIANGLES  = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_REMARKS    = REMARK_LIST');
fprintf(fid1,'%s\n','POINT_KEY_COMPRESSED = NO_LIST');
fprintf(fid1,'%s\n','POINT_KEY_INDICES    = NO_LIST');

fprintf(fid1,'%s%d\n','POINT_NR_LOCATIONS   =  ',num_locs);

fprintf(fid1,'%s\n','POINT_NR_TIMEPTS     =  1');
fprintf(fid1,'%s\n','POINT_TYPE           =  1');
fprintf(fid1,'%s\n','POINT_COORD_SYSTEM   =  0');
fprintf(fid1,'%s\n','POINT_PLOT_FLAGS     =  1');
fprintf(fid1,'%s\n','POINT_PLOT_FLAGS_EX  =  0');
fprintf(fid1,'%s\n','POINT_PLOT_COLOR_1   =  2');
fprintf(fid1,'%s\n','POINT_PLOT_COLOR_2   =  0');
fprintf(fid1,'%s\n','POINT_PLOT_SHAPE     =  4');
fprintf(fid1,'%s\n','POINT_PLOT_SURFACE   =  4');
fprintf(fid1,'%s\n','POINT_PLOT_TRANSPA   =  100');
fprintf(fid1,'%s\n','POINT_PLOT_CLIPPING  =  0');
fprintf(fid1,'%s\n','POINT_PLOT_TEXTSIZE  =  10');
fprintf(fid1,'%s\n','POINT_PLOT_BORDER    =  50');
fprintf(fid1,'%s\n','POINT_PLOT_ADJACENT  =  0');
fprintf(fid1,'%s\n','POINT_PLOT_CLOSED    =  0');
fprintf(fid1,'%s\n','POINT_PLOT_TYPE_1    =  0');
fprintf(fid1,'%s\n','POINT_PLOT_TYPE_2    =  0');
fprintf(fid1,'%s\n','POINT_PLOT_TYPE_3    =  0');
fprintf(fid1,'%s\n','POINT_T_FIRST        =  0');
fprintf(fid1,'%s\n','POINT_T_DELTA        =  0');
fprintf(fid1,'%s\n','POINT_DISTANCE       =  0');
fprintf(fid1,'%s\n','POINT_AREA           =  0');
fprintf(fid1,'%s\n','POINT_VOLUME         =  0');
fprintf(fid1,'%s\n','POINT_SYMBOLSIZE     =  3');
fprintf(fid1,'%s\n','POINT_SYMBOLSCALE    =  3');
fprintf(fid1,'%s\n','POINT_LINEWIDTH      =  1');
fprintf(fid1,'%s\n','POINT_PLOT_DIST_1    =  0');
fprintf(fid1,'%s\n','POINT_PLOT_DIST_2    =  0');
fprintf(fid1,'%s\n','POINT_PLOT_DIST_3    =  0');
fprintf(fid1,'%s\n','POINT_PLOT_PLANE_1   = (0,0,1)');
fprintf(fid1,'%s\n','POINT_PLOT_PLANE_2   = (0,0,1)');
fprintf(fid1,'%s\n','POINT_PLOT_PLANE_3   = (0,0,1)');
fprintf(fid1,'%s\n','POINT_KEYWORDS END');
fprintf(fid1,'%s\n','');
fprintf(fid1,'%s\n','POINT_DESCRIPTION START_LIST	# Do not edit!');
fprintf(fid1,'%s\n','Localize');
fprintf(fid1,'%s\n','(no description available)');
fprintf(fid1,'%s\n','POINT_DESCRIPTION END_LIST');
fprintf(fid1,'%s\n','');
fprintf(fid1,'%s\n','POINT_TRAFO START_LIST	# Do not edit!');
fprintf(fid1,'%s\n','1		 0		 0		 0');
fprintf(fid1,'%s\n','0		-1		 0		 0');
fprintf(fid1,'%s\n','0		 0		-1		 0');
fprintf(fid1,'%s\n','0		 0		 0		 1');
fprintf(fid1,'%s\n','POINT_TRAFO END_LIST');
fprintf(fid1,'%s\n','');
fprintf(fid1,'%s\n','LOCATION_LIST START	# Do not edit!');
fprintf(fid1,'%s\n','LIST_DESCRIPTION     = Locations');
fprintf(fid1,'%s\n','LIST_UNITS           = mm');

fprintf(fid1,'%s%d\n','LIST_NR_ROWS         =  ',num_locs);

fprintf(fid1,'%s\n','LIST_NR_COLUMNS      =  3');
fprintf(fid1,'%s\n','LIST_NR_TIMEPTS      =  1');
fprintf(fid1,'%s\n','LIST_VALID           =  1');
fprintf(fid1,'%s\n','LIST_BINARY          =  0');
fprintf(fid1,'%s\n','LIST_TYPE            =  1');
fprintf(fid1,'%s\n','LIST_TRAFO_TYPE      =  1');
fprintf(fid1,'%s\n','LIST_FIRST_COLUMN    =  1');
fprintf(fid1,'%s\n','LIST_INDEX_MIN       = -1');
fprintf(fid1,'%s\n','LIST_INDEX_MAX       = -1');
fprintf(fid1,'%s\n','LIST_INDEX_ABS_MAX   = -1');
fprintf(fid1,'%s\n','LOCATION_LIST END');
fprintf(fid1,'%s\n','');
fprintf(fid1,'%s\n','LOCATION_LIST START_LIST	# Do not edit!');

for i=1:num_locs
    fprintf(fid1,'%f\t %f\t %f\n',unrot_newtargetL(i,1), unrot_newtargetL(i,2), unrot_newtargetL(i,3));
end

fprintf(fid1,'%s\n','LOCATION_LIST END_LIST');
fprintf(fid1,'%s\n','');
fprintf(fid1,'%s\n','REMARK_LIST START	# Do not edit!');
fprintf(fid1,'%s\n','   LIST_DESCRIPTION     = Remarks');
fprintf(fid1,'%s\n','   LIST_UNITS           = nn');
fprintf(fid1,'%s%d\n','   LIST_NR_ROWS         =  ',num_locs);
fprintf(fid1,'%s\n','   LIST_NR_COLUMNS      =  40');
fprintf(fid1,'%s\n','   LIST_NR_TIMEPTS      =  1');
fprintf(fid1,'%s\n','   LIST_VALID           =  1');
fprintf(fid1,'%s\n','   LIST_BINARY          =  0');
fprintf(fid1,'%s\n','   LIST_TYPE            =  5');
fprintf(fid1,'%s\n','   LIST_TRAFO_TYPE      =  0');
fprintf(fid1,'%s\n','   LIST_FIRST_COLUMN    =  1');
fprintf(fid1,'%s\n','   LIST_INDEX_MIN       = -1');
fprintf(fid1,'%s\n','   LIST_INDEX_MAX       = -1');
fprintf(fid1,'%s\n','   LIST_INDEX_ABS_MAX   = -1');
fprintf(fid1,'%s\n','REMARK_LIST END');
fprintf(fid1,'%s\n','');
fprintf(fid1,'%s\n','REMARK_LIST START_LIST	# Do not edit!');
for i=1:num_locs
    fprintf(fid1,'%s%d\n','Entry ',i);
end
fprintf(fid1,'%s\n','REMARK_LIST END_LIST');
fclose(fid1);









