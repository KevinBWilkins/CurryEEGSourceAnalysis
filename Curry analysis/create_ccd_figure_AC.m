% create_ccd_figure_AC()
% function create_ccd_figure_AC()
% By Albert Chen, Jun Yao
%
% This function read the cortical current density (CCD) information
% in the target local areas and then finds the size of the active area
% in each region of interest. It can use either ANOVA or  
% thresholding to determine the active area.
%
% Last Modified 1/17/07 - AC
%
% Manual changes must be made in the following sections:
%   1. Enter subject parameters before running program, specify when 0ms occurs
%
%
% Sample usage:
% 
% Curry_analysis_AC()
% [Subject] = Curry_analysis_AC(Subject)


function [Subject] = create_ccd_figure_AC(varargin)

if ~isempty(varargin) > 0
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
warning off MATLAB:griddata:DuplicateDataPoints
warning off MATLAB:MKDIR:DirectoryExists

leftshift = 0;
rightshift = 0;
alphashift = 0;

%1. Enter subject parameters before running the program

%*************************************************
%AMPUTEE SUBJECTS
%-------------------------------------------------------------------------
% subject = {'AKpre1'}
% subjects = subject;
% 
% %taskList={'abd_un,'ee_un','hcl_un'};
% taskList={'hcl_un'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 41];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'AKpre1cortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%--------------------------------------------------------------------------
% subject = {'AKpre2'}
% subjects = subject;
% 
% %taskList={'abd','ee','hcl'};
% taskList={'hcl'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 41];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'AKpre2cortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%-------------------------------------------------------------------------
% subject = {'LS1'}
% subjects = subject;
% 
% %taskList={'abd','hcl'};
% taskList={'hcl'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'LS1cortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%-------------------------------------------------------------------------
% subject = {'LS2'}
% subjects = subject;
% 
% %taskList={'ee','ho'};
% taskList={'ho'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'LS2cortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%-------------------------------------------------------------------------
% subject = {'JS2'}
% subjects = subject;
% 
% %taskList={'hcl','ef','ws};
% taskList={'ef'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'JS2cortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% leftshift = 2;
% rightshift = 0;
% alphashift = -1.5;
%-------------------------------------------------------------------------
% subject = {'JS2new'}
% subjects = subject;
% 
% %taskList={'hcl'};
% taskList={'hcl'};
% 
% timePoints_total = [155];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-300.1 302]};     %time ranges of all tasks
% timeRange_BL_total = {[-702.8 -546.4]};   %time ranges of baselines
% plotstop=78;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% 
% im_cortex_name = 'JS2newcortex.jpg';
% 
% leftshift = 2;
% rightshift = 0;
% alphashift = -1.5;
%-------------------------------------------------------------------------
% subject = {'JS3'}
% subjects = subject;
% 
% %taskList={'abd_un','abd'};
% taskList={'abd_un'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'JS3cortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% leftshift = 2;
% rightshift = 0;
% alphashift = -1.5;
%-------------------------------------------------------------------------
% subject = {'JS4'}
% subjects = subject;
% 
% %taskList={'hcl_un'};
% taskList={'hcl_un'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'JS4cortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% leftshift = 2;
% rightshift = 0;
% alphashift = -1.5;
%-------------------------------------------------------------------------
% subject = {'CM1new2'}
% subjects = subject;
% 
% %taskList={'abd'};
% taskList={'abd'};
% 
% timePoints_total = [155];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-300.1 302]};     %time ranges of all tasks
% timeRange_BL_total = {[-702.8 -546.4]};   %time ranges of baselines
% plotstop=78;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% 
% im_cortex_name = 'CM1new2cortex.jpg';
% 
% leftshift = -5;
% rightshift = 3;
%-------------------------------------------------------------------------
% subject = {'CM2new2'}
% subjects = subject;
% 
% %taskList={'hc'};
% taskList={'hc'};
% 
% timePoints_total = [155];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-300.1 302]};     %time ranges of all tasks
% timeRange_BL_total = {[-702.8 -546.4]};   %time ranges of baselines
% plotstop=78;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% 
% im_cortex_name = 'CM2new2cortex.jpg';
% 
% leftshift = -5;
% rightshift = 3;
%-------------------------------------------------------------------------
% subject = {'CM103006'}
% subjects = subject;
% 
% %taskList={'abd_un','hcl_un'};
% 
% taskList={'abd_un'};
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% 
% % taskList={'hcl_un'};
% % timePoints_total = [155];          %number of time points in each task
% % timePoints_BL_total = [41];       %number of time points in each baseline
% % timeRange_total = {[-300.1 302]};     %time ranges of all tasks
% % timeRange_BL_total = {[-702.8 -546.4]};   %time ranges of baselines
% % plotstop=78;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% 
% im_cortex_name = 'CM103006cortex.jpg';
% 
% leftshift = -5;
% rightshift = 3;
%-------------------------------------------------------------------------
% subject = {'CM110206'}
% subjects = subject;
% 
% %taskList={'abd','hcl'};
% taskList={'abd'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'CM110206cortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% leftshift = -5;
% rightshift = 3;


%*************************************************************************
%STROKE SUBJECTS
%-------------------------------------------------------------------------
% subject = {'FNrobot'}
% subjects = subject;
% 
% %taskList={'RE','0RE'};
% taskList={'0RE'};
% 
% timePoints_total = [207];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% 
% im_cortex_name = 'FNrobotcortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%-------------------------------------------------------------------------
% subject = {'FNrobot2'}
% subjects = subject;
% 
% %taskList={'25RE'};
% taskList={'25RE'};
% 
% timePoints_total = [207];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% 
% im_cortex_name = 'FNrobot2cortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%-------------------------------------------------------------------------
% subject = {'WSrobot'}
% subjects = subject;
% 
% %taskList={'0RE'};
% taskList={'0RE'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'WSrobotcortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% leftshift = -4;
%-------------------------------------------------------------------------
% subject = {'WSrobot2'}
% subjects = subject;
% 
% %taskList={'RE'};
% taskList={'RE'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'WSrobot2cortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% leftshift = -4;
%-------------------------------------------------------------------------
% subject = {'WSrobot3'}
% subjects = subject;
% 
% %taskList={'25RE'};
% taskList={'25RE'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'WSrobot3cortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% leftshift = -4;
%-------------------------------------------------------------------------
% subject = {'NTrobot1'}
% subjects = subject;
% 
% % %taskList={'RE','0RE'};
% % taskList={'0RE'};
% % 
% % timePoints_total = [207, 207];          %number of time points in each task
% % timePoints_BL_total = [41, 52];       %number of time points in each baseline
% % timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% % timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% % plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% 
% %taskList={'RE_new'};
% taskList={'RE_new'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-502.5 303.0]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% plotstop=129;
% 
% im_cortex_name = 'NTrobot1cortex.jpg';
% 
% leftshift = 2;
%-------------------------------------------------------------------------
% subject = {'NTrobot2'}
% subjects = subject;
% 
% % %taskList={'25RE'};
% % taskList={'25RE'};
% % 
% % timePoints_total = [207, 207];          %number of time points in each task
% % timePoints_BL_total = [41, 52];       %number of time points in each baseline
% % timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% % timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% % plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% 
% %taskList={'25RE_new'};
% taskList={'25RE_new'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-502.5 303.0]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% plotstop=129;
% 
% im_cortex_name = 'NTrobot2cortex.jpg';
% leftshift = 2;
%-------------------------------------------------------------------------
% subject = {'JWrobot'}
% subjects = subject;
% 
% %taskList={'RE'};
% taskList={'RE'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'JWrobotcortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%-------------------------------------------------------------------------
% subject = {'JWrobot2'}
% subjects = subject;
% 
% %taskList={'0RE','25RE'};
% taskList={'25RE'};
% 
% timePoints_total = [207, 207];          %number of time points in each task
% timePoints_BL_total = [41, 52];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6],[-700 100]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7],[-1950 -1800]};   %time ranges of baselines
% 
% im_cortex_name = 'JWrobot2cortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT

%*************************************************************************
%CONTROL SUBJECTS
%-------------------------------------------------------------------------
% subject = {'YZrobot1'}
% subjects = subject;
% 
% %taskList={'RE_new','0RE_new'};
% taskList={'RE_new'};
% 
% timePoints_total = [207];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-502.5 303.0]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% plotstop=129;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% 
% im_cortex_name = 'YZrobot1cortex.jpg';
%-------------------------------------------------------------------------
% subject = {'YZrobot2'}
% subjects = subject;
% 
% %taskList={'25RE'};
% taskList={'25RE'};
% 
% timePoints_total = [207];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-502.5 303.0]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% plotstop=129;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% 
% im_cortex_name = 'YZrobot2cortex.jpg';
%-------------------------------------------------------------------------
% subject = {'SNrobot'}
% subjects = subject;
% 
% %taskList={'RE','0RE','25RE'};
% taskList={'25RE'};
% 
% timePoints_total = [207];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% 
% im_cortex_name = 'SNrobotcortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%-------------------------------------------------------------------------
% subject = {'GMrobot'}
% subjects = subject;
% 
% % %taskList={'RE','0RE','25RE'};
% % taskList={'25RE'};
% % 
% % timePoints_total = [207];          %number of time points in each task
% % timePoints_BL_total = [41];       %number of time points in each baseline
% % timeRange_total = {[-701.9 103.6]};     %time ranges of all tasks
% % timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% % plotstop=181; %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% 
% %taskList={'RE_new'};
% taskList={'RE_new'};
% 
% timePoints_total = [207];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-502.5 303.0]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% plotstop=129; %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% 
% im_cortex_name = 'GMrobotcortex.jpg'; 
%-------------------------------------------------------------------------
% subject = {'BSrobot1'}
% subjects = subject;
% 
% % %taskList={'RE','0RE','25RE'};
% % taskList={'RE'};
% % 
% % timePoints_total = [207];          %number of time points in each task
% % timePoints_BL_total = [41];       %number of time points in each baseline
% % timeRange_total = {[-701.9 103.6]};     %time ranges of all tasks
% % timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% % plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
% 
% %taskList={'RE_new'};
% taskList={'RE_new'};
% 
% timePoints_total = [207];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-502.5 303.0]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% plotstop=129;
% 
% im_cortex_name = 'BSrobot1cortex.jpg';

%-------------------------------------------------------------------------
% subject = {'JDrobot'}
% subjects = subject;
% 
% %taskList={'RE','0RE','25RE'};
% taskList={'25RE'};
% 
% timePoints_total = [207];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-701.9 103.6]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% 
% im_cortex_name = 'JDrobotcortex.jpg';
% plotstop=181;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT

%*************************************************************************
%Jun's subjects
subjectNames_Jun = ['DW','SN','WS','MLM'];
%-------------------------------------------------------------------------
% subject = {'DW'}
% subjects = subject;
% 
% %taskList={'abd','ef'};
% taskList={'abd'};
% 
% timePoints_total = [155];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-701.9 -99.7]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% 
% im_cortex_name = 'DWcortex.jpg';
% plotstop=155;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%-------------------------------------------------------------------------
% subject = {'SN'}
% subjects = subject;
% 
% %taskList={'abd','ef'};
% taskList={'abd'};
% 
% timePoints_total = [155];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-701.9 -99.7]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% 
% im_cortex_name = 'SNcortex.jpg';
% plotstop=155;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%-------------------------------------------------------------------------
% subject = {'WS'}
% subjects = subject;
% 
% %taskList={'abd','ef'};
% taskList={'ef'};
% 
% timePoints_total = [155];          %number of time points in each task
% timePoints_BL_total = [41];       %number of time points in each baseline
% timeRange_total = {[-701.9 -99.7]};     %time ranges of all tasks
% timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines
% 
% im_cortex_name = 'WScortex.jpg';
% plotstop=155;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT
%-------------------------------------------------------------------------
subject = {'MLM'}
subjects = subject;

%taskList={'abd','ef'};
taskList={'ef'};

timePoints_total = [155];          %number of time points in each task
timePoints_BL_total = [41];       %number of time points in each baseline
timeRange_total = {[-701.9 -99.7]};     %time ranges of all tasks
timeRange_BL_total = {[-1953.1 -1796.7]};   %time ranges of baselines

im_cortex_name = 'MLMcortex.jpg';
plotstop=155;   %where 0ms is       <-----------MUST CHANGE FOR EACH SUBJECT

%2. Some initialization

%individual time plots, what to plot in "for(timestart:timeinc:timestop)"
timestart=1;
timestop=max(timePoints_total);
timeinc=ceil((timestop-timestart+1)/16);  %want about 16 individual plots of 50ms

if strfind(subjectNames_Jun,subject{1})
    timeinc = ceil(207/16);
end

%transformation matrix to coregister curry reconstruction file with cortex image
%COORDINATES ARE ROTATED 18degrees about the x-axis
alpha = (-21+alphashift)*pi/180; 
beta = 0*pi/180;    
gamma = 0*pi/180; 

Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];    %rotation about x axis
Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];        %rotation about y axis
Rz = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];    %rotation about z axis
Rx = Rx*Ry*Rz;

%average cortex plot, what to average in plot in "for plotstart:plotstop"
%usually from -700ms to 0ms or something
%so have to figure out which indices correspond to those times
plotstart=plotstop-timeinc;

subNo=length(subjects);
plotCort=0;
%**************************************************************************

%3. Run through loop of subjects

for k=1:subNo
    subjectName=subject{k};
    mkdir('C:\Albert Chen\matlab\Curry analysis\figures\',subjectName)
    
    %--------------------------------------------------------------------
    %4. Import areas of interest and central sulcus
    
    %names of subjects with special processing procedure (no points selected by Carolina)
    subjectNames_allregion = ['FNrobot', 'FNrobot2', 'YZrobot1', 'YZrobot2', 'LS1', 'LS2','JS2','JS2new','JS3','JS4','GMrobot','CM103006','CM110206','CM1new2','CM2new2'];
    subjectNames_newregion = [];
    %subjectNames_newregion = ['WSrobot', 'WSrobot2', 'WSrobot3'];
    
    if strfind(subjectNames_allregion,subjectName)
        %using Curry 5.0 Albert's allregion pts
        edge1_file_name=[Base_dir,subjectName,'\Points\','allregions_ac.pom'];%allregions
        [edge1,Ecount,ENR]=read_Curry_file4_AC(edge1_file_name,'LOCATION',0,0);
        censul_file_name1=[Base_dir,subjectName,'\Points\','CS_ac.pom'];   %central sulcus
    else
        %using Curry 5.0 Carolina's pts
        edge1_file_name=[Base_dir,subjectName,'\Points\','M1_rt_cc.pom'];%M1_rt
        edge2_file_name=[Base_dir,subjectName,'\Points\','M1_lt_cc.pom'];%M1_lt
        edge3_file_name=[Base_dir,subjectName,'\Points\','S1_rt_cc.pom'];%S1_rt
        edge4_file_name=[Base_dir,subjectName,'\Points\','S1_lt_cc.pom'];%S1_lt
        edge5_file_name=[Base_dir,subjectName,'\Points\','PM_rt_cc.pom'];%PM_rt
        edge6_file_name=[Base_dir,subjectName,'\Points\','PM_lt_cc.pom'];%PM_lt
        edge7_file_name=[Base_dir,subjectName,'\Points\','SMA_rt_cc.pom'];%SMA_rt
        edge8_file_name=[Base_dir,subjectName,'\Points\','SMA_lt_cc.pom'];%SMA_lt
        if strfind(subjectNames_newregion,subjectName)
            edge9_file_name=[Base_dir,subjectName,'\Points\','new_region.pom'];%newregion
            [edge9,Ecount,ENR]=read_Curry_file4_AC(edge9_file_name,'LOCATION',0,0);  
        end
        
        %READ IN LOCATIONS OF POINTS
        [edge1,Ecount,ENR]=read_Curry_file4_AC(edge1_file_name,'LOCATION',0,0);
        [edge2,Ecount,ENR]=read_Curry_file4_AC(edge2_file_name,'LOCATION',0,0);
        [edge3,Ecount,ENR]=read_Curry_file4_AC(edge3_file_name,'LOCATION',0,0);
        [edge4,Ecount,ENR]=read_Curry_file4_AC(edge4_file_name,'LOCATION',0,0);
        [edge5,Ecount,ENR]=read_Curry_file4_AC(edge5_file_name,'LOCATION',0,0);
        [edge6,Ecount,ENR]=read_Curry_file4_AC(edge6_file_name,'LOCATION',0,0);
        [edge7,Ecount,ENR]=read_Curry_file4_AC(edge7_file_name,'LOCATION',0,0);
        [edge8,Ecount,ENR]=read_Curry_file4_AC(edge8_file_name,'LOCATION',0,0);
        
        %Read in central sulcus file, may be different name  
        censul_file_name1=[Base_dir,subjectName,'\Points\','CS_cc.pom'];   %central sulcus
    end
    
    [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name1,'LOCATION',0,0);
    censul2 = censul1;
    censul=censul1;
    
    %----------------------------------------------------------------------
    %5. Get indices of mask-selected points over all areas of interest
    
    [ROI_ind,ROI_ind_group,Rx,allmasks]=select_ROI_AC(subjectName,-21);
    
    %----------------------------------------------
    %6. Combine points to make areas such as SMA, PM, M1, and S1
    
    %     areanames = {'SMA','PM','M1','S1','ContraHem','IpsiHem','Whole'};
    areanames = {'SMA_rt','PM_rt','M1_rt','S1_rt','ContraHem','IpsiHem','Whole'};
    %     areanames = {'SMA_lt','PM_lt','M1_lt','S1_lt','ContraHem','IpsiHem','Whole'};
    
    if strfind(subjectNames_allregion,subjectName)
        areas{7} = edge1; %when using allregions
    else
        areas{1} = [edge7; edge8]; %SMA
        areas{2} = [edge5; edge6]; %PM
        areas{3} = [edge1; edge2]; %M1
        areas{4} = [edge3; edge4]; %S1
        
        %     areas{1} = [edge8];   %Ipsilateral areas only
        %     areas{2} = [edge6];
        %     areas{3} = [edge2];
        %     areas{4} = [edge4];
        
        %     areas{1} = [edge7];     %Contralateral areas only
        %     areas{2} = [edge5];
        %     areas{3} = [edge1];
        %     areas{4} = [edge3];
        
        areas{5} = [edge1;edge3;edge5;edge7];       %RIGHT HEMISPHERE
        areas{6} = [edge2;edge4;edge6;edge8];       %LEFT HEMISPHERE
        areas{7} = [edge1;edge3;edge5;edge7;edge2;edge4;edge6;edge8];   %ALL MOTOR AREAS
        
        if strfind(subjectNames_newregion,subjectName)
            areas{7} = [edge1;edge2;edge3;edge4;edge5;edge6;edge7;edge8;edge9];   %ALL AREAS
        end
    end
    %end of picking areas of interest
    %---------------------------------------------------------------

    %7. Loop through tasks- such as RE, 0RE, 25RE, etc
    for i=1:length(taskList)
        cur_task=taskList{i};
        
        timePoints = timePoints_total(i);
        timePoints_BL = timePoints_BL_total(i);
        timeRange = timeRange_total{i};
        timeRange_BL = timeRange_BL_total{i};
        
        %---------------------------------------------------------------
        %8. Loop through areas of interest- 
        
        %1=M1, 2=S1, 3=PM, 4=SMA, 5=right hem, 6=left hem, 7=whole cortex
        for areanum = 7:7
            fprintf(1, 'processing areanum %d... \n',areanum)
            edge = areas{areanum};
            
            %----------------------------------
            cdr_file=[cur_task,'.cdr'];
            cdr_BL_file=[cur_task,'_BL.cdr'];
            
            disp('********************');
            disp(cdr_file)
            cdr_file_name=[Base_dir,subjectName,'\new\',cdr_file];
            cdr_BL_file_name=[Base_dir,subjectName,'\new\',cdr_BL_file];%BASELINE
            
            %--------------------------------------------
            %9. Get strengths of points in target areas
            %
            %9a. get all cortex locations from *.cdr file
            fprintf(1, 'get cortex locations... ')
            [cortexL,Lcount,LNR]=read_Curry_file3(cdr_file_name,'LOCATION',0,0);
            fprintf(1, 'done\n')
            
            %9b. get all cortex strengths corresponding closest to edge points
            fprintf(1, 'get region strengths... ')
            [cortexC,targetL,targetC]=find_cortex_strengths_Jun(cdr_file_name,cortexL,ROI_ind,timePoints,0);
            fprintf(1, 'done\n')
            
            if ~isempty(ROI_ind_group)
                fprintf(1, 'get left_region strengths... ')
                [cortexC,left_targetL,left_targetC]=find_cortex_strengths_Jun(cdr_file_name,cortexL,ROI_ind(ROI_ind_group>0),timePoints,0);
                fprintf(1, 'done\n')
                
                %figure out which coordinates of targetL belong to right and left
                %zeros refer to right ROI, 1's refer to left ROI
                targetL_group = ismember(targetL,left_targetL);
                
                left_group = targetL_group>0;
                left_group = left_group(:,1);
                left_group = left_group>0;
                right_group = ~left_group;
                
                fprintf(1, 'get right_region strengths... ')
                right_targetC = targetC(right_group,:);
                fprintf(1, 'done\n')
            end
            
            %9c. get baseline strengths of points in target areas
            fprintf(1, 'get cortex baseline strengths... ')
            if ~strfind(subjectNames_Jun,subjectName)
                [cortexC_BL,targetL_BL,targetC_BL]=find_cortex_strengths_Jun(cdr_BL_file_name,cortexL,ROI_ind,timePoints_BL,0);
            end
                fprintf(1, 'done\n')
            
            if ~isempty(ROI_ind_group) & ~strfind(subjectNames_Jun,subjectName)
                fprintf(1, 'get right_region strengths... ')
                right_targetC_BL = targetC_BL(right_group,:);
                fprintf(1, 'done\n')
                
                fprintf(1, 'get left_region strengths... ')
                left_targetC_BL = targetC_BL(left_group,:);
                fprintf(1, 'done\n')
            end
            %--------------------------------------------------------------
            
            %10. Rotate cortex and target points
            cortexL=(Rx*cortexL')';
            targetL=(Rx*targetL')';
            
            if ~isempty(ROI_ind_group)
                rot_right_targetL=targetL(right_group,:);
                rot_left_targetL=targetL(left_group,:);
                
                figure(2)
                hold on;
                plot3(rot_right_targetL(:,1),rot_right_targetL(:,2),rot_right_targetL(:,3),'r.', 'MarkerSize',6);
                plot3(rot_left_targetL(:,1),rot_left_targetL(:,2),rot_left_targetL(:,3),'g.', 'MarkerSize',6);  
            end
            %--------------------------------------------------------------
            
            if ~strfind(subjectNames_Jun,subjectName)
                %11. Baseline and cortex strength amplitude processing

                %11a. Baseline processing- Jun's method
                targetC_BL1=mean(targetC_BL,2);

                %targetC takes last number of points equal to winPoints
                %             targetC=targetC(:,(end-winPoints+1):end);
                winPoints = plotstop-plotstart+1;
                targetC1=targetC(:,plotstart:plotstop);

                %subtract baseline
                targetC1=targetC1-repmat(targetC_BL1,1,winPoints);
                %             targetC1=targetC1-repmat(targetC_BL1,1,plotstop-plotstart+1);
                %             targetC1=targetC1(:,plotstart:plotstop);
                
            else    %don't use baseline
                winPoints = plotstop-plotstart+1;
                targetC1=targetC(:,plotstart:plotstop);
            end
            
            [locationNo,TimeNo]=size(targetC1);

            targetC1_sum=sum(targetC1,2);
            targetC1_norm=targetC1_sum/max(targetC1_sum);

            %--------Jun added on 09/25/06 for the case that some
            %targetC_BL is larger than the targetC
            targetC1_norm(find(targetC1_norm<0))=0;

            meantargetC1_norm = mean(targetC1_norm)
            [muhat, muci] = expfit(targetC1_norm)
            ind=find(targetC1_norm>muci(2));

            ActiveRatio1(k,i,areanum)=length(ind)/locationNo; %active area ratio

            figure(70+i)
            subplot(1,2,1)
            e_y = sort(targetC1_norm);
            plot(e_y,'Color','g','Linewidth',2);
            hold on;
            e_x = (1:length(targetC1_norm))';
            line([1 e_x(end)],[muci(1) muci(1)],'Color','r','Linewidth',2);
            line([1 e_x(end)],[muci(2) muci(2)],'Color','r','Linewidth',2);


            log_e_y = log(e_y);
            %fit = polyfit(e_x(500:1300), log_e_y(500:1300), 1);
            fit = polyfit(e_x, log_e_y, 1);

            %             plot(e_x, exp(fit(2)).*exp(fit(1)*e_x),'Color','b','Linewidth',2);

            %             figure(90+i);plot(log_e_y)
            %             hold on;
            %             line([400 400],[-10 max(log_e_y)],'Color','r','Linewidth',2);
            %             line([1200 1200],[-10 max(log_e_y)],'Color','r','Linewidth',2);

            
            %--------------------------------------------------------------
            %11b. Baseline processing- Albert's method
            
            targetC2 = targetC(:,1:plotstop);
            
            if ~strfind(subjectNames_Jun,subjectName)
                totaltime_BL = timeRange_BL(2)-timeRange_BL(1);
                times1_BL = timeRange_BL(1):totaltime_BL/(timePoints_BL-1):timeRange_BL(2);

                meantargetC_BL = mean(targetC_BL,2);
                sortmeantargetC_BL = sort(meantargetC_BL);

                %Subtract baseline strengths from all strengths
                targetC2 = targetC2 - repmat(meantargetC_BL,1,plotstop-plotstart+1);
            end
                
            targetC2(find(targetC2<0))=0;
            
            maxtargetC = max(max(targetC2(:,1:plotstop)));
           
            [m,n] = size(targetC2);
            sorttargetC = sort(reshape(targetC2,m*n,1))/maxtargetC;
            
            %Figure out threshold
            %[muhat, muci] = expfit(sortmeantargetC_BL);
            meansorttargetC = mean(sorttargetC)
            [muhat2, muci2] = expfit(sorttargetC)
            ind=find(sorttargetC>muci2(2));

            %Threshold = max(sortmeantargetC_BL)
            Threshold = muci2(2);
            
            %thresholds = [.85; .80; .75; .70; .65; .60; .50; .40; .20]*maxtargetC;
            
            PercentThresholdofMax = Threshold
            
            %Plot of thresholds
            figure(70+i)
            subplot(1,2,2)
            %plot(sortmeantargetC_BL,'Color','m','Linewidth',2)
            hold on;
            plot(sorttargetC,'Color','g','Linewidth',2)
            line([1 length(sorttargetC)],[Threshold Threshold],'Color','r','Linewidth',2)
            %line([1 length(sorttargetC)],[thresholds thresholds],'Color','m','Linewidth',2)
            line([1 length(sorttargetC)],[muci2(1) muci2(1)],'Color','r','Linewidth',2);
            line([1 length(sorttargetC)],[muci2(2) muci2(2)],'Color','r','Linewidth',2);
            
            axis([1 length(sorttargetC) 0 1])
            
            ind=find(targetC1_norm>muci(2));
            ActiveRatio1(k,i,areanum)=length(ind)/locationNo; %active area ratio
            
            thresholds = Threshold;
            
            
            %------------------------------------------
            %12. Plot pretty pictures
            samp_inc = 0.8;
            
            %ROTATED CORTEX LOCATIONS TO MATCH TOP DOWN LOOK
            rot_cortexL = cortexL;
            mincortX=min(rot_cortexL(:,1));
            maxcortX=max(rot_cortexL(:,1));
            mincortY=min(rot_cortexL(:,2));
            maxcortY=max(rot_cortexL(:,2));
            bigx=mincortX:samp_inc:maxcortX;
            bigy=mincortY:samp_inc:maxcortY;
            [bigXI,bigYI] = meshgrid(bigx,bigy);
            l_bigy = length(bigy);      %dimensions of big cortex locations
            l_bigx = length(bigx);
            
            %Rotated target locations
            %find limits of rotated target locations
            rot_targetL = targetL;
            minrot_X=min(rot_targetL(:,1));
            maxrot_X=max(rot_targetL(:,1));
            minrot_Y=min(rot_targetL(:,2));
            maxrot_Y=max(rot_targetL(:,2));
            minrot_Z=min(rot_targetL(:,3));
            maxrot_Z=max(rot_targetL(:,3));
            
            %change these limits to be same as coordinates from cortex
            [y_dum,i_dum] = min(abs(bigx-minrot_X));
            minrot_X = bigx(i_dum);
            [y_dum,i_dum] = min(abs(bigx-maxrot_X));
            maxrot_X = bigx(i_dum);
            [y_dum,i_dum] = min(abs(bigy-minrot_Y));
            minrot_Y = bigy(i_dum);
            [y_dum,i_dum] = min(abs(bigy-maxrot_Y));
            maxrot_Y = bigy(i_dum);
            
            rot_x=minrot_X:samp_inc:maxrot_X;
            rot_y=minrot_Y:samp_inc:maxrot_Y;
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
            rot_edgex=minrot_edgeX:samp_inc:maxrot_edgeX;
            rot_edgey=minrot_edgeY:samp_inc:maxrot_edgeY;
            [rot_edgeXI,rot_edgeYI] = meshgrid(rot_edgex,rot_edgey); 
            rot_edgeZI = griddata(rot_edgeL(:,1),rot_edgeL(:,2),rot_edgeL(:,3),rot_edgeXI,rot_edgeYI);
            
                       
            %13. Get highest z-values at each x-y coordinate (only top surface of cortex)
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
            
            if ~strfind(subjectNames_Jun,subjectName)
                top_cortexC_BL = cortexC_BL(indzs,:);
            end
            
%             figure
%             plot3(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),'LineWidth',7,'color','r','LineStyle','.')
%             hold on;
%             plot3(top_rot_cortexL(:,1),top_rot_cortexL(:,2),top_rot_cortexL(:,3),'LineWidth',7,'color','b','LineStyle','.')
%             plot3(rot_targetL(:,1),rot_targetL(:,2),rot_targetL(:,3),'LineWidth',7,'color','g','LineStyle','.')
%             view(0,0);
%             bigZI = griddata(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),bigXI,bigYI);
%             bigZI = griddata(top_rot_cortexL(:,1),top_rot_cortexL(:,2),top_rot_cortexL(:,3),bigXI,bigYI);
            
            %ROTATE everything else too
            rot_censul = (Rx*censul')';
            
            %------------------------------------------
            %13b. Find grid locations for sampled points on total grid
            %and create a mask where 1 = region, 0 = outside region
            
            %determine whether mask has already been created
%             if length(varargin) < 1
%                 createmaskflag = 1;
%             else
%                 size_allmasks = size(allmasks);
%                 if size_allmasks(2) < areanum
%                     createmaskflag = 1;
%                 else
%                     if isempty(allmasks{i,areanum})
%                         createmaskflag = 1;
%                     else
%                         createmaskflag = 0;
%                     end
%                 end
%             end
            if isempty(allmasks)
                createmaskflag = 1;
            else
                createmaskflag = 0;
            end
            
            
            %14. IMAGE PROCESSING OF CORTEX PICTURE
            %used to overlay region of interest onto picture of cortex
            im_cortex_file_name=[Base_dir,subjectName,'\Points\',im_cortex_name];%image file location
            imA = imread(im_cortex_file_name);       %raw picture
            imB = im2bw(imA,.95);               %make black and white
            se = strel('square',5);             
            imC = imclose(imB,se);              %close up holes
            imD = ~imC;                         %flip white and black
            imE = imfill(imD,'holes');          %fill holes
            [im_i,im_j] = find(imE>0);          %find filled area
            imF = imE(min(im_i):max(im_i),min(im_j):max(im_j));     %
            imG = bwperim(imF,8);
            
            imA_resize = imA((min(im_i)-3):(max(im_i)+1),(min(im_j)-12+leftshift):(max(im_j)+16-rightshift),:);
%             imA_resize = imA((min(im_i)-3):(max(im_i)+1),(min(im_j)-16):(max(im_j)+17),:);
%             imA_resize = imA(min(im_i):max(im_i),min(im_j):max(im_j),:);
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
            
            %get actual perimeter locations around cortex
            imB2 = im2bw(imA_resize2,.95);               %make black and white
            imC2 = imclose(imB2,se);              %close up holes
            imD2 = ~imC2;                         %flip white and black
            imE2 = imfill(imD2,'holes');          %fill holes
            imG2 = bwperim(imE2,8);
            [imG_pts_i,imG_pts_j] = find(imG2>0);
            for ind_cort=1:5:length(imG_pts_i)
                loc_cortex2(ind_cort,:) = [x_locs(imG_pts_j(ind_cort)),y_locs(imG_pts_i(ind_cort))];
            end
            loc_cortex2 = loc_cortex2(1:5:end,:);
            
            
            %determine midline
            %midline1_x = x_locs(round(length(x_locs)/2));
            [P,S] = polyfit(loc_cortex(:,2),loc_cortex(:,1),1);
            midline1_x = P(1)*y_locs+P(2);
            %[P,S] = polyfit(rot_cortexL(:,2),rot_cortexL(:,1),1)
            %midline2_x = P(1)*y_locs+P(2);
            
            %Figure displaying cortex picture, border of cortex, central sulcus, midline                
            figure(1)
            imshow(imA_resize2,'XData',[mincortX maxcortX],'YData',[maxcortY mincortY])
            axis on
            axis([mincortX maxcortX mincortY maxcortY])
            axis xy
            hold on;
            scatter3(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),7,'y');
            scatter3(loc_cortex2(:,1),loc_cortex2(:,2),zeros(length(loc_cortex2(:,1)),1),7,'m')
            plot3(rot_censul(:,1),rot_censul(:,2),rot_censul(:,3),'Markersize',15,'color','b','LineStyle','.')
            %scatter3(rot_censul(:,1),rot_censul(:,2),rot_censul(:,3),3,'g');
            %scatter3(rot_targetL(:,1),rot_targetL(:,2),rot_targetL(:,3),7,'r');
            
            %plot midline
            %line([midline1_x midline1_x],[mincortY maxcortY],'Color','red','Linestyle','--')
%             plot(midline1_x,y_locs,'color','m','Linewidth',2,'Linestyle','--')
            %plot(midline2_x,y_locs,'color','m','Linewidth',2,'Linestyle','--')
            
            view(0,90)
            axis on
            axis equal
            axis([mincortX maxcortX mincortY maxcortY])
            
            %resize imA_resize2 so that pixel size = samp_inc
            imA_resize3 = imresize(imA_resize2,[length(bigy),length(bigx)],'bilinear');
            
            
            createmaskflag = 0;
            %15. CREATE MASK OF ROTATED TARGET LOCATIONS IF NO MASK
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
                allmasks{i,areanum} = flipud(bigmask4);
                changemaskflag = 1;
                %close figure
                close(100)
            else
                changemaskflag = 0;
                %DISPLAY SAVED MASK
                figure(100)
                h1 = imagesc(imA_resize3);
                hold on;
                
                bigmask4 = flipud(allmasks);
                h2 = imagesc(bigmask4,[0 3]);
                colormap vga
                set(h2,'AlphaData',.5);

                but = 1;
                changemaskflag = 0;
                mask_value = 1;
                while but == 1 | but == 3 | but == 2 | but == 49 | but == 50
                    [xi,yi,but] = ginput(1);
                    if but==1 | but == 3 | but == 2 | but == 49 | but == 50
                        if but == 49
                            mask_value = 1;
                        elseif but == 50
                            mask_value = 2;
                        end

                        if round(yi) < (l_bigy-2) && round(yi) > 2 && round(xi) < (l_bigx-2) && round(xi) > 2
                            if but == 1
                                bigmask4((round(yi)-2):(round(yi)+2),(round(xi)-2):(round(xi)+2)) = mask_value;
                                %also change mask4
                            elseif but == 2
                                bigmask4((round(yi)-2):(round(yi)+2),(round(xi)-2):(round(xi)+2)) = 2;
                            elseif but == 3
                                bigmask4((round(yi)-2):(round(yi)+2),(round(xi)-2):(round(xi)+2)) = 0;
                                %also change mask4
                            end
                            %imagesc(mask4+2*edge_mask);
                            set(h2,'CData',bigmask4)
                        end
                    end
                end
                allmasks = flipud(bigmask4);

                close(100)
%                 changemaskflag = 0;
            end
            
            fprintf(1, 'done\n')       
            
            bigmask4 = allmasks;
            mask_perim = bwperim(bigmask4,8);
            mask_perim_nums = union(mask_perim_nums,[areanum]);
            mask_perims{areanum} = mask_perim;
            
            %-----------------------------------
            
            
            %-------------------------------------
            %Make new edge and find new targetC and targetL locations
            
            fprintf(1, 'create new mask... ')
            
            %----------------------------------
            %15. CONVERT MASKS to LOCATION MASKS
            
            %make location mask for ROI
            [mask_m,mask_n] = find(bigmask4>0);
            Xmat = repmat(rot_x,length(rot_y),1);
            Ymat = repmat(rot_y,length(rot_x),1)';
            
            for mask_ind = 1:length(mask_m)
                
                %if inside target locations area
                if bigx(mask_n(mask_ind))<maxrot_X & bigx(mask_n(mask_ind))>minrot_X & bigy(mask_m(mask_ind))<maxrot_Y & bigy(mask_m(mask_ind))>minrot_Y
                    %figure out where mask index would be
                    
                    dx = bigx(mask_n(mask_ind)); 
                    dy = bigy(mask_m(mask_ind));
                    [minval,minrow] = min((dx-Xmat).^2 + (dy-Ymat).^2);
                    [minval,mincol] = min(minval);
                    
                    loc_y = minrow(mincol); 
                    loc_x = mincol;
                    dx2 = rot_x(loc_x); 
                    dy2 = rot_y(loc_y);
                    
                    rot_new_edge(mask_ind,:) = [dx2 dy2 rot_ZI(loc_y,loc_x)];
                else
                    %if outside target locations area
                    rot_new_edge(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
                end
            end
            
            %------------------------------------
            %make location mask for the perimeter of the ROI
            [mask_m,mask_n] = find(mask_perim>0);
            
            for mask_ind = 1:length(mask_m)
                
                %if inside target locations area
                if bigx(mask_n(mask_ind))<maxrot_X && bigx(mask_n(mask_ind))>minrot_X && bigy(mask_m(mask_ind))<maxrot_Y && bigy(mask_m(mask_ind))>minrot_Y
                    %figure out where mask index would be
                    
                    dx = bigx(mask_n(mask_ind)); 
                    dy = bigy(mask_m(mask_ind));
                    [minval,minrow] = min((dx-Xmat).^2 + (dy-Ymat).^2);
                    [minval,mincol] = min(minval);
                    
                    loc_y = minrow(mincol); 
                    loc_x = mincol;
                    dx2 = rot_x(loc_x); 
                    dy2 = rot_y(loc_y);
                    
                    rot_new_perim(mask_ind,:) = [dx2 dy2 rot_ZI(loc_y,loc_x)];
                else
                    %if outside target locations area
                    rot_new_perim(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
                end
            end
            
            regions{i,areanum} = rot_new_perim;
            alledges{i,areanum} = rot_new_edge;
            fprintf(1, 'done\n')
            %--------------------------------------------------------------
            
            %16. GET NEW STRENGTHS AND BASELINE STRENGTHS if mask created or changed
            if changemaskflag
                
                if ~isempty(ROI_ind_group)
                    fprintf(1, 'get new strengths... ')
                    
                    [left_newtargetL,left_newtargetC,new_left_group]=find_toplayer_cortex_strengths_AC(top_rot_cortexL,top_cortexC,rot_left_targetL,rot_edgeL);
                    [right_newtargetL,right_newtargetC,new_right_group]=find_toplayer_cortex_strengths_AC(top_rot_cortexL,top_cortexC,rot_right_targetL,rot_edgeL);
                    fprintf(1, 'done\n')
                    
                    fprintf(1, 'get new baseline strengths...')
                    %[left_newtargetL_BL,left_newtargetC_BL,Tind]=find_toplayer_cortex_strengths_AC(top_rot_cortexL,top_cortexC_BL,rot_left_targetL,rot_edgeL);
                    %[right_newtargetL_BL,right_newtargetC_BL,Tind]=find_toplayer_cortex_strengths_AC(top_rot_cortexL,top_cortexC_BL,rot_right_targetL,rot_edgeL);
                    fprintf(1, 'done\n')
                    
                    left_group = new_left_group;
                    right_group = new_right_group;
                    newtargetL = [left_newtargetL;right_newtargetL];
                    newtargetL_BL = [left_newtargetL_BL;right_newtargetL_BL];
                    newtargetC = [left_newtargetC;right_newtargetC];
                    newtargetC_BL = [left_newtargetC_BL;right_newtargetC_BL];
                        
                else
                    fprintf(1, 'get new strengths... ')
                    %[newtargetL,newtargetC]=find_new_cortex_strengths_AC(top_rot_cortexL,top_cortexC,rot_new_edge,timePoints,0);
                    [newtargetL,newtargetC,Tind]=find_toplayer_cortex_strengths_AC(top_rot_cortexL,top_cortexC,rot_new_edge,rot_edgeL);
                    fprintf(1, 'done\n')
                    
                    fprintf(1, 'get new baseline strengths...')
                    %[newtargetL_BL,newtargetC_BL,Tind]=find_toplayer_cortex_strengths_AC(top_rot_cortexL,top_cortexC_BL,rot_new_edge,rot_edgeL);
                    fprintf(1, 'done\n')
                end
                
            else
                newtargetL = targetL;
                newtargetC = targetC;
                %newtargetL_BL = targetL_BL;
                %newtargetC_BL = targetC_BL;
            end
            
            %---------------------------------
            %16b. New baseline processing- Albert's method
            
            newtargetC2 = newtargetC;
            
            if ~strfind(subjectNames_Jun,subjectName)
                meantargetC_BL = mean(newtargetC_BL,2);
                sortmeantargetC_BL = sort(meantargetC_BL);
                meantargetC_BL2 = mean(newtargetC(:,plotstart:(plotstart+4*timeinc)),2);

                %Subtract baseline strengths from all strengths
                newtargetC2 = newtargetC2 - repmat(meantargetC_BL,1,size(newtargetC2,2));
            end
            newtargetC2(find(newtargetC2<0))=0;
            
            %subtract first 50 ms of CDR baseline activity from all strengths
%             newtargetC2 = newtargetC2 - repmat(meantargetC_BL2,1,size(newtargetC2,2));
%             newtargetC2(find(newtargetC2<0))=0;
            
            
            [m,n] = size(newtargetC2);
            sorttargetC = sort(reshape(newtargetC2,m*n,1));
            
            %take top 5% of sources as maxtargetC
            numtop5 = round(length(sorttargetC)/20);
            
            %maxtargetC = max(max(newtargetC2(:,plotstart:plotstop)));
            realmaxtargetC = max(max(newtargetC2(:,1:plotstop)))
            maxtargetC = mean(sorttargetC(end-numtop5:end));
            maxvalues(areanum) = maxtargetC
            
            %normalize sorttargetC
            sorttargetC = sorttargetC/maxtargetC;
            
            %Figure out threshold
            %[muhat, muci] = expfit(sortmeantargetC_BL);
            meansorttargetC = mean(sorttargetC)
            [muhat2, muci2] = expfit(sorttargetC)
            ind=find(sorttargetC>muci2(2));

            newThreshold = muci2(2);
            newThresholds(areanum) = newThreshold;
            %thresholds = [.85; .80; .75; .70; .65; .60; .50; .40; .20]*maxtargetC;
            thresholds = newThreshold;
            
            PercentThresholdofMax = newThreshold
            percentmaxthresholds(areanum) = PercentThresholdofMax;
            
            figure(71+i)
            %plot(sortmeantargetC_BL,'Color','m','Linewidth',2)
            hold on;
            plot(sorttargetC,'Color','g','Linewidth',2)
            line([1 length(sorttargetC)],[newThreshold newThreshold],'Color','r','Linewidth',2)
            %line([1 length(sorttargetC)],[muci2(1) muci2(1)],'Color','r','Linewidth',2);
            %line([1 length(sorttargetC)],[muci2(2) muci2(2)],'Color','r','Linewidth',2);
            
            axis([1 length(sorttargetC) 0 1])                    
            
            %normalize targetC
            %newtargetC = newtargetC2/maxtargetC;
            newtargetC = newtargetC2;
            
            %-----------------------------------------------------------
            %17. Plot one whole averaged cortex over time period chosen
            figure(22+4*areanum+1)
            hold on
            timenow=plotstart:plotstop;
            C = mean(newtargetC(:,timenow),2);
            maxmeanC = max(max(C));
            fprintf('max mean strength = %d\n',maxmeanC)
            
            
            %---------------------------------------------------------
            %unrotate points back, so that can create a .pom file with new locations
%             unRx = [1 0 0; 0 cos(-alpha) sin(-alpha); 0 -sin(-alpha) cos(-alpha)];
%             unRy = [cos(-beta) 0 -sin(-beta); 0 1 0; sin(-beta) 0 cos(-beta)];
%             unRz = [cos(-gamma) sin(-gamma) 0; -sin(-gamma) cos(-gamma) 0; 0 0 1];
%             unRx = unRx*unRy*unRz;
%             unrot_newtargetL = (unRx*newtargetL')';
            
            %---------------------------------------
            %create new .pom file of new locations
%             write_new_pom_file(areanames,areanum,Base_dir,subjectName,unrot_newtargetL);
            
            
            %------------------------------------
            %find limits of new target locations
            newminrot_X=min(newtargetL(:,1));
            newmaxrot_X=max(newtargetL(:,1));
            newminrot_Y=min(newtargetL(:,2));
            newmaxrot_Y=max(newtargetL(:,2));
            newminrot_Z=min(newtargetL(:,3));
            newmaxrot_Z=max(newtargetL(:,3));
            
            %3/15/07- need to take subset of rot_x and rot_y, not new coordinates
            %with same sampling frequency
            [y_dum,i_dum] = min(abs(rot_x-newminrot_X));
            newminrot_X = rot_x(i_dum);
            [y_dum,i_dum] = min(abs(rot_x-newmaxrot_X));
            newmaxrot_X = rot_x(i_dum);
            [y_dum,i_dum] = min(abs(rot_y-newminrot_Y));
            newminrot_Y = rot_y(i_dum);
            [y_dum,i_dum] = min(abs(rot_y-newmaxrot_Y));
            newmaxrot_Y = rot_y(i_dum);    
            
            newrot_x=[newminrot_X:samp_inc:newmaxrot_X];
            newrot_y=[newminrot_Y:samp_inc:newmaxrot_Y];
            
            [newrot_XI,newrot_YI] = meshgrid(newrot_x,newrot_y); 
            newrot_ZI = griddata(newtargetL(:,1),newtargetL(:,2),newtargetL(:,3),newrot_XI,newrot_YI);
            
            CI = griddata(newtargetL(:,1),newtargetL(:,2),C,newrot_XI,newrot_YI);
            
            %SURFACE PLOT OVERLAYED WITH MASK PERIM
            surface('XData',newrot_XI,'YData',newrot_YI,'ZData',newrot_ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6)
            
            plot3(rot_new_perim(:,1),rot_new_perim(:,2),rot_new_perim(:,3),'LineWidth',7,'color','r','LineStyle','.')
            
            %scatter3(rot_censul(:,1),rot_censul(:,2),rot_censul(:,3),3,'g');
            
            view(0,90)
            axis off
            axis equal
            colorbar
            v_caxis = caxis;
            %close(22+4*areanum+1)
            
            %---------------           
            %need mask for CI
            [mp_i,mp_j] = find(mask_perim>0);
            new_mask_perim = mask_perim(min(mp_i):max(mp_i),min(mp_j):max(mp_j));
            new_mask_perim_filled = imfill(new_mask_perim,'holes');
            %new_mask_perim_filled = imfill(mask_perim,'holes');     
            
            newminrot_X=min(newrot_x);  %use these as mins and maxes because samp_inc may be too large
            newmaxrot_X=max(newrot_x);
            newminrot_Y=min(newrot_y);
            newmaxrot_Y=max(newrot_y);
            
            [numrow_CI,numcol_CI] = size(CI);
            [i_CI,j_CI] = find(CI>0);
            %need to figure out how many rows and columns of NaN's there are on borders of CI
            %min(i_CI) tells you where first row non-NaN entry of CI is
            %max(i_CI) tells you where last row non-NaN entry of CI is
            %min(j_CI) tells you where first column non-NaN entry of CI is
            %max(j_CI) tells you where last column non-NaN entry of CI is
            
            %so pad new_mask_perim_filled with zeros on all sides
            
            %pad min(j_CI)-1 columns in front
            %pad numcol_CI-max(j_CI) columns to end
            %pad min(i_CI)-1 rows to front
            %pad numrow_CI-max(i_CI) rows to end
            
            %match filled perimeter mask locations to CI locations to make mask for CI
            new_mask_perim_filled_resize = new_mask_perim_filled;
            if (min(rot_new_perim(:,1)) < newrot_x(min(j_CI)) && min(rot_new_perim(:,1)) < newminrot_X)    %if perim mask is to left of first CI NaN and actual
                disp 'a'
                %just take mask to the right of newminrot_X
                new_mask_perim_filled_resize = new_mask_perim_filled_resize(:,(1+round((newminrot_X-min(rot_new_perim(:,1)))/samp_inc)):(end));
            elseif min(rot_new_perim(:,1)) > newrot_x(min(j_CI)) || (min(rot_new_perim(:,1)) <= newrot_x(min(j_CI)) && min(rot_new_perim(:,1)) > newminrot_X)    %if perim mask is to right of first CI actual or is inbetween NaN and actual
                disp 'b'
                %need to pad first n columns with zero, where n = diff between newminrot_X and min_x of mask
                nmpfr_size = size(new_mask_perim_filled_resize);
                new_mask_perim_filled_resize_dum = zeros(nmpfr_size(1),nmpfr_size(2)+(round((min(rot_new_perim(:,1))-newminrot_X)/samp_inc)));
                new_mask_perim_filled_resize_dum(:,(1+(round((min(rot_new_perim(:,1))-newminrot_X)/samp_inc))):(end)) = new_mask_perim_filled_resize;
                new_mask_perim_filled_resize = new_mask_perim_filled_resize_dum;
            end
            
            if max(rot_new_perim(:,1)) < newrot_x(max(j_CI)) || (max(rot_new_perim(:,1)) >= newrot_x(max(j_CI)) && max(rot_new_perim(:,1)) < newmaxrot_X)    %if perim mask is to left of last CI actual or is inbetween actual and NaN
                disp 'c'
                %pad zeros to right of max_x of rot_new_perim
                new_mask_perim_filled_resize(:,(end+1):(end+round((newmaxrot_X-max(rot_new_perim(:,1)))/samp_inc))) = zeros(size(new_mask_perim_filled,1),round((newmaxrot_X-max(rot_new_perim(:,1)))/samp_inc));
            elseif max(rot_new_perim(:,1)) > newmaxrot_X   %if perim mask is to right of last CI actual and NaN
                disp 'd'
                %just take mask to the left of newmaxrot_X
                new_mask_perim_filled_resize = new_mask_perim_filled_resize(:,1:(end-round((max(rot_new_perim(:,1))-newmaxrot_X)/samp_inc)));
            end
            %-------------
            
            if min(rot_new_perim(:,2)) < newminrot_Y   %if perim mask is below lowest CI row actual and NaN
                disp 'e'
                %take rows above newminrot_y
                new_mask_perim_filled_resize = new_mask_perim_filled_resize((1+round((newminrot_Y-min(rot_new_perim(:,2)))/samp_inc)):(end),:);
            elseif min(rot_new_perim(:,2)) > newrot_y(min(i_CI)) || (min(rot_new_perim(:,2)) <= newrot_y(min(i_CI)) && min(rot_new_perim(:,2)) > newminrot_Y)
                %pad first few rows with zeros, then shift everything up to top
                disp 'f'
                new_mask_perim_filled_resize_new = zeros(size(new_mask_perim_filled_resize,1)+round((min(rot_new_perim(:,2))-newminrot_Y)/samp_inc),size(new_mask_perim_filled_resize,2));
                new_mask_perim_filled_resize_new(1+(round((min(rot_new_perim(:,2))-newminrot_Y)/samp_inc)):end,:) = new_mask_perim_filled_resize;
                new_mask_perim_filled_resize = new_mask_perim_filled_resize_new;
            end
            
            if (max(rot_new_perim(:,2)) < newrot_y(max(i_CI))) || (max(rot_new_perim(:,2)) >= newrot_y(max(i_CI)) && max(rot_new_perim(:,2)) < newmaxrot_Y)    %if perim mask is below highest row of CI actual or is inbetween highest CI actual and NaN
                disp 'g'
                %pad top rows with zeros
                new_mask_perim_filled_resize((end+1):(end+round((newmaxrot_Y-max(rot_new_perim(:,2)))/samp_inc)),1:size(new_mask_perim_filled_resize,2)) = zeros(round((newmaxrot_Y-max(rot_new_perim(:,2)))/samp_inc),size(new_mask_perim_filled_resize,2));
            elseif max(rot_new_perim(:,2)) > newmaxrot_Y     %if perim mask is higher than all rows
                disp 'h'
                %take rows below to newmaxrot_Y
                new_mask_perim_filled_resize = new_mask_perim_filled_resize(1:(end-round((max(rot_new_perim(:,2))-newmaxrot_Y)/samp_inc)),:);
            end
            
%             if min(rot_new_perim(:,1)) < newrot_x(min(j_CI))
%                 disp 'a'
%                 %just take mask to the right of newminrot_X
%                 new_mask_perim_filled_resize = new_mask_perim_filled_resize(:,(1+round((newminrot_X-min(rot_new_perim(:,1)))/samp_inc)):(end));
%             elseif min(rot_new_perim(:,1)) > newrot_x(min(j_CI))
%                 disp 'b'
%                 %need to pad first n rows with zero, where n = diff between newminrot_X and min_x of mask
%                 nmpfr_size = size(new_mask_perim_filled_resize);
%                 new_mask_perim_filled_resize_dum = zeros(nmpfr_size(1),nmpfr_size(2)+(round((min(rot_new_perim(:,1))-newminrot_X)/samp_inc)));
%                 new_mask_perim_filled_resize_dum(:,(1+(round((min(rot_new_perim(:,1))-newminrot_X)/samp_inc))):(end)) = new_mask_perim_filled_resize;
%                 new_mask_perim_filled_resize = new_mask_perim_filled_resize_dum;
%             end
%             
%             if max(rot_new_perim(:,1)) < newrot_x(max(j_CI))
%                 disp 'c'
%                 %pad zeros to right of max_x of rot_new_perim
%                 new_mask_perim_filled_resize(:,(end+1):(end+round((newmaxrot_X-max(rot_new_perim(:,1)))/samp_inc))) = zeros(size(new_mask_perim_filled,1),round((newmaxrot_X-max(rot_new_perim(:,1)))/samp_inc));
%             elseif max(rot_new_perim(:,1)) > newrot_x(max(j_CI))
%                 disp 'd'
%                 %just take mask to the left of newmaxrot_X
%                 new_mask_perim_filled_resize = new_mask_perim_filled_resize(:,1:(end-round((max(rot_new_perim(:,1))-newmaxrot_X)/samp_inc)));
%             end
%             
%             if min(rot_new_perim(:,2)) < newrot_y(min(i_CI))
%                 disp 'e'
%                 %take rows above newminrot_y
%                 new_mask_perim_filled_resize = new_mask_perim_filled_resize((1+round((newminrot_Y-min(rot_new_perim(:,2)))/samp_inc)):(end),:);
%             elseif min(rot_new_perim(:,2)) > newrot_y(min(i_CI))
%                 %pad first few rows with zeros, then shift everything up to
%                 %top
%                 disp 'f'
%                 new_mask_perim_filled_resize_new = zeros(size(new_mask_perim_filled_resize,1)+round((min(rot_new_perim(:,2))-newminrot_Y)/samp_inc),size(new_mask_perim_filled_resize,2));
%                 new_mask_perim_filled_resize_new((round((min(rot_new_perim(:,2))-newminrot_Y)/samp_inc)):(end-1),:) = new_mask_perim_filled_resize;
%                 new_mask_perim_filled_resize = new_mask_perim_filled_resize_new;
%             end
%             
%             if max(rot_new_perim(:,2)) < newrot_y(max(i_CI))
%                 disp 'g'
%                 %pad top rows with zeros
%                 new_mask_perim_filled_resize((end+1):(end+round((newmaxrot_Y-max(rot_new_perim(:,2)))/samp_inc)),1:size(new_mask_perim_filled_resize,2)) = zeros(round((newmaxrot_Y-max(rot_new_perim(:,2)))/samp_inc),size(new_mask_perim_filled_resize,2));
%             elseif max(rot_new_perim(:,2)) > newrot_y(max(i_CI))
%                 disp 'h'
%                 %take rows below newmaxrot_Y
%                 new_mask_perim_filled_resize = new_mask_perim_filled_resize(1:(end-round((max(rot_new_perim(:,2))-newmaxrot_Y)/samp_inc)),:);
%             end

            nmpfr_size = size(new_mask_perim_filled_resize);
            if ~isequal(nmpfr_size,size(CI))
                disp 'new_mask_perim_filled_resize is not same size as CI'
                %pause
            end
            new_mask_perim_filled_resize = double(new_mask_perim_filled_resize);
            new_mask_perim_filled_resize(find(new_mask_perim_filled_resize>0)) = .6;
            
             
            %-------------------------
            
            targetCI = newtargetC;
            
            if ~isempty(ROI_ind_group)
                if (ROI_ind_group(1) == 1)
                    left_targetCI = targetCI(1:sum(left_group>0),:);
                    right_targetCI = targetCI((sum(left_group>0)+1):end,:);
                    left_newtargetL = newtargetL(1:sum(left_group>0),:);
                    right_newtargetL = newtargetL((sum(left_group>0)+1):end,:);
                else
                    right_targetCI = targetCI(1:sum(right_group>0),:);
                    left_targetCI = targetCI((sum(right_group>0)+1):end,:);
                    right_newtargetL = newtargetL(1:sum(right_group>0),:);
                    left_newtargetL = newtargetL((sum(right_group>0)+1):end,:);
                end
            end
            
            figure(50)
            scatter3(left_newtargetL(:,1),left_newtargetL(:,2),left_newtargetL(:,3),'g');
            hold on;
            scatter3(right_newtargetL(:,1),right_newtargetL(:,2),right_newtargetL(:,3),'m');
            
            %---------------------------------------------------
            %18a. Figure out statistical threshold of activation
            
            %                 [m,n] = size(targetCI);
            %                 longtargetCI = reshape(targetCI,m*n,1);
            %                 sorttargetCI = sort(longtargetCI);
            %                 
            %                 V1='B'; %baseline
            %                 V2='E'; %activation (execution)
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
            %18b. Statistical determination of significant differences
            %from baseline over range of time window
            
            %                 %number of points
            %                 %[locationNo_BL,TimeNo_BL]=size(targetC_BL);
            %                 [locationNo,TimeNo]=size(targetC);
            %                 
            %                 %loop through all points and determine whether significant
            %                 % difference from baseline over the range of time window
            %                 siglocationInd = [];
            %                 ActCount = 0;
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
            
            %----------------------------------------------------------
            %19. CCDR quantification
            
            fprintf(1, 'quantifying results... ')
            
            totaltime = timeRange(2)-timeRange(1);
            times1 = timeRange(1):totaltime/(timePoints-1):timeRange(1)+50;
            
            %Figure out locations of active areas within time bins
            
            figure(50+4*areanum)
            set(50+4*areanum,'Color',[1 1 1])
            index=0;
%             
%             %8 subplots
%             for timenow = timestart+8*timeinc:timeinc:timestop
            %16 subplots
            for timenow = timestart:timeinc:timestop
                index=index+1;
                
                %8 or 16 subplots
                subplot(4,4,index)
                %subplot(2,4,index)
                imshow(imA_resize2,'XData',[mincortX maxcortX],'YData',[maxcortY mincortY])
                axis([mincortX maxcortX mincortY maxcortY])
                axis off
                hold on
                
                timestart_seg = timenow;
                timestop_seg = timenow + timeinc-1;
                if timestop_seg > timestop
                    timestop_seg = timestop;
                end
                
                %SURFACE PLOT of active areas, all same amplitude
                C = mean(newtargetC(:,timestart_seg:timestop_seg),2);
%                 [thr_i,thr_j] = find(C<newThreshold);
%                 Thr_ind = sub2ind(size(C),thr_i,thr_j);
%                 C_bin = ones(size(C));
%                 C_bin(Thr_ind) = 0;
                C_bin = double(C>newThreshold);
                CI_bin = griddata(newtargetL(:,1),newtargetL(:,2),C_bin,newrot_XI,newrot_YI);
                h3 = surface('XData',newrot_XI,'YData',newrot_YI,'ZData',newrot_ZI,'CData',CI_bin,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6);
                set(h3,'FaceAlpha','flat','AlphaData',new_mask_perim_filled_resize);
                set(gca,'ALim',[0 1]);
                
                view(0,90)
                TR = timeRange_total{i};
                title([num2str(TR(1)+(TR(2)-TR(1))*(timestart_seg-1)/(timePoints_total(i)-1),5),' to ',num2str(TR(1)+(TR(2)-TR(1))*(timestop_seg-1)/(timePoints_total(i)-1),5),'ms'],'Fontsize',14,'FontWeight','bold')
                axis off
                caxis([0 1])
                axis equal
                axis xy
                axis tight
            end
            hgsave(gcf,[pwd,'\figures\',subjectName,'\',subjectName,'_',cur_task,'_active.fig']);
            
            %Quantify active area vs time
            %number of locations in that area above threshold value
            overthresholds{i,areanum} = sum(targetCI>newThreshold,1);
            
            %percentage of active area that is over threshold
            totalarea = size(targetCI,1);
            
            percentactiveareas{i,areanum} = overthresholds{i,areanum}/totalarea;
            
            %total sum of strengths of active areas over threshold at each time point
            totalstrengths{i,areanum} = sum(targetCI.*double(targetCI>newThreshold),1);
            
            %average strength of active areas over threshold at each time point
            %dummy_avgstrengths = totalstrengths{oind,i,areanum,methodInd}./overthresholds{oind,i,areanum,methodInd}./maxvalue{i,areanum,methodInd};
            dummy_avgstrengths = totalstrengths{i,areanum}./overthresholds{i,areanum};
            dummy_avgstrengths(isnan(dummy_avgstrengths)) = 0;
            avgstrengths{i,areanum} = dummy_avgstrengths;
            
            if ~isempty(ROI_ind_group)
                left_overthresholds{i,areanum} = sum(left_targetCI>newThreshold,1);
                right_overthresholds{i,areanum} = sum(right_targetCI>newThreshold,1);
                left_totalarea = size(left_targetCI,1);     
                right_totalarea = size(right_targetCI,1); 
                left_percentactiveareas{i,areanum} = left_overthresholds{i,areanum}/left_totalarea;
                right_percentactiveareas{i,areanum} = right_overthresholds{i,areanum}/right_totalarea;
                left_totalstrengths{i,areanum} = sum(left_targetCI.*double(left_targetCI>newThreshold),1);
                right_totalstrengths{i,areanum} = sum(right_targetCI.*double(right_targetCI>newThreshold),1);
                
                %hemisphere laterality index is weighted measure of strength and location
                wHLIleft = left_totalstrengths{i,areanum}.*left_overthresholds{i,areanum};
                wHLIright = right_totalstrengths{i,areanum}.*right_overthresholds{i,areanum};
                wHLIs{i,areanum} = (wHLIright-wHLIleft)./(wHLIright+wHLIleft);
                %unweighted hemisphere laterality index
                HLIleft = left_overthresholds{i,areanum};   %unweighted HLI
                HLIright = right_overthresholds{i,areanum}; %unweighted HLI
                HLIs{i,areanum} = (HLIright-HLIleft)./(HLIright+HLIleft);
                
            end
            %"center of activity" is strength weighted locations, averaged together 
            %strengths = overthresholds{oind,i,areanum,methodInd}
            %locations = targetL;
            [m_target,n_target] = size(targetCI);
            
%             for col = 1:n_target
%                 S_col = targetCI(:,col);
%                 sum_S_col = sum(S_col,1);
%                 S_col_rep = repmat(S_col,1,3);
%                 wt_loc = newtargetL.*S_col_rep;
%                 mean_loc(:,col) = (1/sum_S_col)*sum(wt_loc,1)';
%             end
            
            %only get locations that are active
            for coord = 1:3
                targetL_M = repmat(newtargetL(:,coord),1,n_target);
                weighted_locs = targetL_M.*double(targetCI>newThreshold).*targetCI;
                mean_locs(coord,1:n_target) = sum(weighted_locs,1)./totalstrengths{i,areanum};
                %replace NaN's with mean coordinate
                mean_locs(isnan(mean_locs)) = mean(newtargetL(:,coord));                
            end
            
            single_mean_loc = mean(mean_locs,2)';
            Subject.single_mean_loc = single_mean_loc;
            Subject.mean_locs = mean_locs;
            
            %plot where center of activity occurs
            figure(101)
            hold on;
            if isequal(cur_task,'RE') | isequal(cur_task,'RE_new')
                m_color = 'b';
                imshow(imA_resize2,'XData',[mincortX maxcortX],'YData',[maxcortY mincortY])
            elseif isequal(cur_task,'0RE') | isequal(cur_task,'0RE_new')
                m_color = 'g';
                imshow(imA_resize2,'XData',[mincortX maxcortX],'YData',[maxcortY mincortY])
            elseif isequal(cur_task,'25RE') | isequal(cur_task,'25RE_new')
                m_color = 'm';
                imshow(imA_resize2,'XData',[mincortX maxcortX],'YData',[maxcortY mincortY])
            else
                m_color = 'm';
                imshow(imA_resize2,'XData',[mincortX maxcortX],'YData',[maxcortY mincortY])
            end
            axis on
            axis([mincortX maxcortX mincortY maxcortY])
            axis xy
            axis equal
            hold on;
            
            %plot3(rot_new_perim(:,1),rot_new_perim(:,2),rot_new_perim(:,3),'LineWidth',7,'color','r','LineStyle','.')
            plot(rot_new_perim(:,1),rot_new_perim(:,2),'LineWidth',7,'color','r','LineStyle','.')
            %scatter3(mean_loc(1,:),mean_loc(2,:),mean_loc(3,:),10,0:2/(length(mean_loc(1,:))-1):2);
            %plot3(single_mean_loc(1),single_mean_loc(2),single_mean_loc(3),'Markersize',30,'color',m_color,'LineStyle','.')
            plot(single_mean_loc(1),single_mean_loc(2),'Markersize',30,'color',m_color,'LineStyle','.')
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
            axis([mincortX maxcortX mincortY maxcortY])
            
            %title(['distance to midline = ',num2str(dist_to_midline)])
            
            %----------------------------------------------------------
            %Display surface on top of image of cortex
            
            
            figure(102)
            set(gcf,'Color',[1 1 1]);
            %top-down image of cortex
            imshow(imA_resize2,'XData',[mincortX maxcortX],'YData',[maxcortY mincortY])
            axis on
            axis([mincortX maxcortX mincortY maxcortY])
            hold on;
            h3 = surface('XData',newrot_XI,'YData',newrot_YI,'ZData',newrot_ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6);
            
            set(h3,'FaceAlpha','flat','AlphaData',new_mask_perim_filled_resize);
            set(gca,'ALim',[0 1]);
            
            %plot new perimeter
%             %new_mask_perim_filled_resize,newrot_XI,newrot_YI
%             se = strel('square',5);             
%             imC3 = imclose(new_mask_perim_filled_resize,se);              %close up holes
%             imE3 = imfill(imC3,'holes');          %fill holes
%             imG3 = bwperim(imE3,8);
%             [imG_pts_i,imG_pts_j] = find(imG3>0);
%             
%             targetL_xs = newrot_XI(1,:);
%             targetL_ys = newrot_YI(:,1);
%             for ind_cort=1:length(imG_pts_i)
%                 rot_new_perim2(ind_cort,:) = [targetL_xs(imG_pts_j(ind_cort)),targetL_ys(imG_pts_i(ind_cort))];
%             end
%             rot_new_perim2 = rot_new_perim2(1:end,:);
%             plot3(rot_new_perim2(:,1),rot_new_perim2(:,2),zeros(length(rot_new_perim2(:,1)),1),'LineWidth',7,'color','r','LineStyle','.');
            %plot3(rot_new_perim(:,1),rot_new_perim(:,2),rot_new_perim(:,3),'LineWidth',7,'color','r','LineStyle','.')
            
            view(0,90)
            axis equal
            axis xy
            caxis(v_caxis);
            axis on;axis equal
            
            a1 = gca;
            a2 = axes;
            caxis(v_caxis);
            colorbar
            set(a2,'Position',get(a1,'Position'));
            axis off;
            
            figure(103)
            set(gcf,'Color',[1 1 1]);
            %top-down image of cortex
            imshow(imA_resize2,'XData',[mincortX maxcortX],'YData',[maxcortY mincortY])
            axis on
            axis([mincortX maxcortX mincortY maxcortY])
            hold on;
            h3 = surface('XData',newrot_XI,'YData',newrot_YI,'ZData',newrot_ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6);
            
            set(h3,'FaceAlpha','flat','AlphaData',new_mask_perim_filled_resize);
            set(gca,'ALim',[0 1]);
            
            %plot new perimeter
%             %new_mask_perim_filled_resize,newrot_XI,newrot_YI
%             se = strel('square',5);             
%             imC3 = imclose(new_mask_perim_filled_resize,se);              %close up holes
%             imE3 = imfill(imC3,'holes');          %fill holes
%             imG3 = bwperim(imE3,8);
%             [imG_pts_i,imG_pts_j] = find(imG3>0);
%             
%             targetL_xs = newrot_XI(1,:);
%             targetL_ys = newrot_YI(:,1);
%             for ind_cort=1:length(imG_pts_i)
%                 rot_new_perim2(ind_cort,:) = [targetL_xs(imG_pts_j(ind_cort)),targetL_ys(imG_pts_i(ind_cort))];
%             end
%             rot_new_perim2 = rot_new_perim2(1:end,:);
%             plot3(rot_new_perim2(:,1),rot_new_perim2(:,2),zeros(length(rot_new_perim2(:,1)),1),'LineWidth',7,'color','r','LineStyle','.');
            plot3(rot_new_perim(:,1),rot_new_perim(:,2),rot_new_perim(:,3),'LineWidth',7,'color','r','LineStyle','.')
            plot3(newtargetL(:,1),newtargetL(:,2),newtargetL(:,3),'LineWidth',7,'color','y','LineStyle','.')
            
            view(0,90)
            axis equal
            axis xy
            caxis(v_caxis);
            axis on;axis equal
            
            a1 = gca;
            a2 = axes;
            caxis(v_caxis);
            colorbar
            set(a2,'Position',get(a1,'Position'));
            axis off;
            %----------------------------------------------------------
            %display normalized image with perim
            
            norm_new_perim(:,1) = (rot_new_perim(:,1)-(mincortX+maxcortX)/2)/((maxcortX-mincortX)/2);
            norm_new_perim(:,2) = (rot_new_perim(:,2)-(mincortY+maxcortY)/2)/((maxcortY-mincortY)/2);
            %norm_new_perim(:,3) = (rot_new_perim(:,3)-(mincortZ+maxcortZ)/2)/(maxcortZ-mincortZ)/2;
            norm_loc_cortex2(:,1) = (loc_cortex2(:,1)-(mincortX+maxcortX)/2)/((maxcortX-mincortX)/2);
            norm_loc_cortex2(:,2) = (loc_cortex2(:,2)-(mincortY+maxcortY)/2)/((maxcortY-mincortY)/2);
            norm_rot_censul(:,1) = (rot_censul(:,1)-(mincortX+maxcortX)/2)/((maxcortX-mincortX)/2);
            norm_rot_censul(:,2) = (rot_censul(:,2)-(mincortY+maxcortY)/2)/((maxcortY-mincortY)/2);
            
            figure(104)
            %top-down image of cortex
            imshow(imA_resize2,'XData',[-1 1],'YData',[1 -1])
            axis on
            axis([-1 1 -1 1])
            hold on;
            
            plot3(norm_new_perim(:,1),norm_new_perim(:,2),rot_new_perim(:,3),'LineWidth',7,'color','r','LineStyle','.')
            scatter3(norm_loc_cortex2(:,1),norm_loc_cortex2(:,2),zeros(length(norm_loc_cortex2(:,1)),1),7,'m')
            plot3(norm_rot_censul(:,1),norm_rot_censul(:,2),rot_censul(:,3),'Markersize',15,'color','b','LineStyle','.')
            
            view(0,90)
            axis equal
            axis xy
            caxis(v_caxis);
            axis on;axis equal
            
            
            
            
            %set NaN entries to 0 if no values above threshold
%             dummy_avg = avgstrengths{i,areanum};
%             dummy_avg(isnan(dummy_avg)) = 0;
            
            %subtract baseline from strengths
            %basestrengths{oind,i,areanum,methodInd} = mean(dummy_avg(1:length(times1)));
            %avgstrengths{oind,i,areanum,methodInd} = (dummy_avg - basestrengths{oind,i,areanum,methodInd});
            
            %maxv = maximum value seen
            maxv = maxtargetC;
            meanSIs{i,areanum} = mean(avgstrengths{i,areanum});
            SIs{i,areanum} = max(avgstrengths{i,areanum});
            
            dummy_avg = avgstrengths{i,areanum};
            
            timestart_seg = plotstop-2*timeinc;
            timestop_seg = plotstop;
            disp([num2str(TR(1)+(TR(2)-TR(1))*(timestart_seg-1)/(timePoints_total(i)-1),5),' to ',num2str(TR(1)+(TR(2)-TR(1))*(timestop_seg-1)/(timePoints_total(i)-1),5),'ms'])
            
            meanSIs_pre = mean(dummy_avg(:,timestart_seg:timestop_seg))
            meanSIs_post = mean(dummy_avg(:,timestop_seg+2*timeinc))
            
            %active area index
            meanAAIs{i,areanum} = mean(percentactiveareas{i,areanum});
            AAIs{i,areanum} = max(percentactiveareas{i,areanum});
            
            dummy_area = percentactiveareas{i,areanum};
            meanAAIs_pre = mean(dummy_area(:,timestart_seg:timestop_seg))
            meanAAIs_post = mean(dummy_area(:,timestop_seg+2*timeinc))
            
            %subtract baseline from areas
            %dummy_area = percentactiveareas{oind,i,areanum,methodInd};
            %baseareas{oind,i,areanum,methodInd} = mean(dummy_area(1:length(times1)));
            %percentactiveareas{oind,i,areanum,methodInd} = (dummy_area - baseareas{oind,i,areanum,methodInd});
            
            if ~isempty(ROI_ind_group)
                dummywHLI = wHLIs{i,areanum};
                meanwHLIs_pre = mean(dummywHLI(:,timestart_seg:timestop_seg))
                meanwHLIs_post = mean(dummywHLI(:,timestop_seg+2*timeinc))
                dummywHLI(isnan(dummywHLI))=0;
                meanwHLIs{i,areanum} = mean(dummywHLI);
                
                dummyHLI = HLIs{i,areanum};
                meanHLIs_pre = mean(dummyHLI(:,timestart_seg:timestop_seg))
                meanHLIs_post = mean(dummyHLI(:,timestop_seg+2*timeinc))
                dummyHLI(isnan(dummyHLI))=0;
                meanHLIs{i,areanum} = mean(dummyHLI);
            end
            fprintf(1, 'done\n')
            
            %20. Plot cortex activity within time bins of 50 ms
            if areanum==7
                
                %20a. Plot                     
                figure(20+4*areanum)
                set(20+4*areanum,'Color',[1 1 1])
                figure(21+4*areanum)
                set(21+4*areanum,'Color',[1 1 1])
                
                index=0;
                %8 subplots
                
                %if starts at -2000, then add 8*timeinc
                %else starts at -1000, add 4*timeinc
                if timePoints_total(1)==155
                    addtime1 = 4*timeinc;
                    addtime2 = -4*timeinc;
                else
                    addtime1 = 8*timeinc;
                    addtime2 = 0;
                end
                
                for timenow = timestart+addtime1:timeinc:timestop+addtime2
                    index=index+1;
                    
                    %8 or 16 subplots
                    figure(20+4*areanum)
                    %subplot(4,4,index)
                    subplot(2,4,index)
                    imshow(imA_resize2,'XData',[mincortX maxcortX],'YData',[maxcortY mincortY])
                    axis([mincortX maxcortX mincortY maxcortY])
                    axis off
                    hold on
                    
                    timestart_seg = timenow;
                    timestop_seg = timenow + timeinc-1;
                    if timestop_seg > timestop
                        timestop_seg = timestop;
                    end
                    
                    mean_loc = mean(mean_locs(:,timestart_seg:timestop_seg),2);
                    norm_mean_locs(1,index) = (mean_loc(1)-(mincortX+maxcortX)/2)/((maxcortX-mincortX)/2);
                    norm_mean_locs(2,index) = (mean_loc(2)-(mincortY+maxcortY)/2)/((maxcortY-mincortY)/2);
                    
                                
                    %SURFACE PLOT
                    C = mean(newtargetC(:,timestart_seg:timestop_seg),2);
                    CI = griddata(newtargetL(:,1),newtargetL(:,2),C,newrot_XI,newrot_YI);
                    h3 = surface('XData',newrot_XI,'YData',newrot_YI,'ZData',newrot_ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6);
                    set(h3,'FaceAlpha','flat','AlphaData',new_mask_perim_filled_resize);
                    set(gca,'ALim',[0 1]);
                    
                    
                    view(0,90)
                    TR = timeRange_total{i};
                    title([num2str(TR(1)+(TR(2)-TR(1))*(timestart_seg-1)/(timePoints_total(i)-1),5),' to ',num2str(TR(1)+(TR(2)-TR(1))*(timestop_seg-1)/(timePoints_total(i)-1),5),'ms'],'Fontsize',14,'FontWeight','bold')
                    axis off
                    caxis([0 1])
                    axis equal
                    axis xy
                    axis tight
                    
                    %plot center of activity
                    
                    figure(21+4*areanum)
                    subplot(2,4,index)
                    
                    hold on;
                    
                    if isequal(cur_task,'RE') | isequal(cur_task,'abd') | isequal(cur_task,'abd_un')
                        m_color = 'b';
                        imshow(imA_resize2,'XData',[-1 1],'YData',[1 -1])
                    elseif isequal(cur_task,'0RE') 
                        m_color = 'g';
                        imshow(imA_resize2,'XData',[-1 1],'YData',[1 -1])
                    elseif isequal(cur_task,'25RE')
                        m_color = 'm';
                        imshow(imA_resize2,'XData',[-1 1],'YData',[1 -1])
                    else
                        m_color = 'm';
                        imshow(imA_resize2,'XData',[-1 1],'YData',[1 -1])
                    end
                    hold on;
                    plot3(norm_mean_locs(1,index),norm_mean_locs(2,index),newmaxrot_Z,'Markersize',20,'color',m_color,'LineStyle','.')
                    %plot(mean_loc(1),mean_loc(2),'Markersize',20,'color','m','LineStyle','.')
                    axis([-1 1 -1 1])

                    title([num2str(TR(1)+(TR(2)-TR(1))*(timestart_seg-1)/(timePoints_total(i)-1),5),' to ',num2str(TR(1)+(TR(2)-TR(1))*(timestop_seg-1)/(timePoints_total(i)-1),5),'ms'],'Fontsize',14,'FontWeight','bold')
                    axis equal
                    axis xy
                    axis tight
                end
                figure(20+4*areanum)
                hgsave(gcf,[pwd,'\figures\',subjectName,'\',subjectName,'_',cur_task,'_8brains.fig']);
            end
                       
             
        end %end of areanum
        %------------------------------------------------------------------
            
    end % end of task
    %----------------------------------------------------------------------
    
end % end of subject Number

beep

%save certain variables and results
Subject.maxvalues = maxvalues;
Subject.newThresholds = newThresholds;
Subject.percentmaxthresholds = percentmaxthresholds;
Subject.avgstrengths = avgstrengths;
Subject.SIs = SIs;
Subject.meanSIs = meanSIs;
Subject.percentactiveareas = percentactiveareas;
Subject.AAIs = AAIs;
Subject.meanAAIs = meanAAIs;

if ~isempty(ROI_ind_group)
    Subject.HLIs = HLIs;
    Subject.meanHLIs = meanHLIs;
end

Subject.allmasks = allmasks;
Subject.mask_perims = mask_perims;
Subject.mask_perim_nums = mask_perim_nums;
Subject.regions = regions;

cog.norm_mean_locs = norm_mean_locs;
cog.loc_cortex2 = loc_cortex2;
cog.imA_resize2 = imA_resize2;
cog.norm_new_perim = norm_new_perim;
cog.norm_rot_censul = norm_rot_censul;


save [pwd,'\cog\',subjectName,'_','cog.mat'] cog 








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

function [targetL,targetC,Tind] = find_toplayer_cortex_strengths_AC(cortexL,contrib,new_edge,old_edge)

fpos=0;

%--------Select the points inside the edge---------------
Tind=[];
[n,tmp]=size(new_edge);
[m,tmp]=size(cortexL);

for i=1:n
    dis2total_edge =  sqrt(sum((old_edge-repmat(new_edge(i,:),size(old_edge,1),1)).^2,2));
    
    if min(dis2total_edge)<7
        dis2cortex= sqrt(sum((cortexL-repmat(new_edge(i,:),size(cortexL,1),1)).^2,2));
        
        [yind,ind] = min(dis2cortex);
        if dis2cortex(ind)<5
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
fid1 = fopen([Base_dir,subjectName,'\Points\','new_',pom_name,'.pom'],'w');

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









