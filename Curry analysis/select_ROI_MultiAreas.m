% function [censul,targetL,targetC_BL,targetC,XI,YI,ZI]=select_ROI (cortexL,targetL,edge)                
function [ROI_ind,ROI_ind_group,ROI_ind_group2,Rx,allmasks]=select_ROI_AC_Jun_MultiAreas(varargin) 
% This function returns the rotated cortex locations and the index of the region of interests. The input
% varargin can be:
%    1. subjectName
%    2. deg_alpha
%    3. remakemask
%    4. SubjectGroup
%    5. AreaNames
% Example:
% subjectName = 'Name';
% deg_alpha = -21;
% remakemask = 1;
% [ROI_ind,ROI_ind_group,Rx,allmasks]=select_ROI_AC_new('DW',-21,[],'control', {'M1_rt_cc';'M1_lt_cc';'S1_rt_cc';'S1_lt_cc';'PM_rt_cc';'PM_lt_cc';'SMA_rt_cc';'SMA_lt_cc';'preSMA_rt_cc';'preSMA_lt_cc';'Pa_rt_cc';'Pa_lt_cc'});
% Please note that: it is important to keep the order always the same for
% all the subjects. (Should take care of this in the code)
warning off MATLAB:griddata:DuplicateDataPoints

subjectName = varargin{1};
deg_alpha = varargin{2};
remakemask = 1;
SubjectGroup = varargin{4};
AreaNames = varargin{5};

timePoints=1;
areanum=1;
yj_debug=1;

%Note: regions are referred to as "edges" in this file
right_edge = [];
left_edge = [];
allmasks = [];
% Base_dir='F:\data\inverse_results\amputee\';
%Base_dir=['F:\data\inverse_results\',SubjectGroup,'\'];
Base_dir=['C:\Users\Kevin\Desktop\LP=1_Results\'];
%if old mask exists in directory
if exist ([Base_dir,subjectName,'\Points\','Mask.mat'],'file')
%     createmaskflag = 1;
    allmasks = [];
    ROI_ind_group = [];
    ROI_ind_group2 = [];
    load ([Base_dir,subjectName,'\Points\','Mask.mat'])
    
    if length(varargin)>2
        remakemask = varargin{3};
    else
        remakemask = 0;
    end
end
    
if remakemask
%     createmaskflag=0;
    

    %----------------------------------------------------------------------
    %Get source locations
    
    %filenames of CDR files
%     cdr_file=['rbic_4'];
%     cdr_file=['abd-full'];
    cdr_file=['LP=1_Open-Average_Pre-Intervention-Open'];

%     cdr_BL_file=['RE_new_BL'];
    
    disp('********************');
    disp(cdr_file)
    cdr_file_name=[Base_dir,subjectName,'\new\',cdr_file,'.cdr'];
    
    
    %EDGE IS USED IN THIS FILE TO MEAN "REGION"
    total_edge = [];    
    for i_area = 1: length (AreaNames)
        %using Curry 5.0 Carolina's pts
        cur_area = AreaNames {i_area};
        cur_edge_file_name = [Base_dir,subjectName,'\Points\',cur_area,'.pom'];
        
        %READ IN LOCATIONS OF POINTS
        if exist (cur_edge_file_name) 
            [edge{i_area},Ecount,ENR]=read_Curry_file4_AC(cur_edge_file_name,'LOCATION',0,0);
        else
            edge{i_area} =[];
        end
        
        total_edge = cat(1,total_edge,edge{i_area});

        %Read in central sulcus file, may be different name
        censul_file_name1=[Base_dir,subjectName,'\Points\','CS-KW.pom'];   %central sulcus
    end
    left_edge = [edge{2};edge{4};edge{6};edge{8}];
    right_edge = [edge{1};edge{3};edge{5};edge{7}];
    areanum = 1;
    %--------------------------------------------------------------
    %if only one central sulcus file
    [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name1,'LOCATION',0,0);
    censul2 = censul1;
    censul=[censul1];
    %----------------------------------------------------------------------    
    %get all cortex locations
    fprintf(1, 'get cortex locations... ')
    [cortexL_mri_voxel,Lcount,LNR,TM]=read_Curry_file3_TM(cdr_file_name,'LOCATION',0,0);
%     [cortexL,cortexL_full,Lcount,LNR]=read_Curry_file3_AC(cdr_file_name,'LOCATION',0,0);
    fprintf(1, 'done\n')

    %if TM is not identity matrix, then transform cortex coordinates
    %because they are for MRI, not Curry
    TM2 = TM(2:4,1);    %translation
    TM1 = TM(2:4,2:4);  %rotation and scaling
    inv_TM1 = inv(TM1);

    %cortexL_mri = (inv_TM1*(cortexL - repmat(TM2',LNR(1),1))')';
    cortexL_curry = (TM1*cortexL_mri_voxel')' + repmat(TM2',LNR(1),1);
    %mirror and rotate 180 degrees- negate x and y coordinates
    cortexL_curry = [-cortexL_curry(:,1:2) cortexL_curry(:,3)];     %this matches curry mri locations to curry real world locations
    cortexL = cortexL_curry;        %replace cortex mri coordinates with curry coordinates
        
    %--------------------------------------------
    %get all cortex locations corresponding closest to edge points
    fprintf(1, 'get ROI locations... ')
    %use matlab picked ROI regions for total_edge if you want to just modify an old mask
    if exist ([Base_dir,subjectName,'\Points\','Mask.mat'],'file')
        
        targetL = find_cortex_locations_AC(cortexL,total_edge);
        
        
    else    %find cortex locations nearest to points picked in Curry
        targetL = find_cortex_locations_AC(cortexL,total_edge);
    end
    
    fprintf(1, 'done\n')

    %--------------------------------------------
    %preparation for mask selection- rotate, choose top layer, image processing, etc.
    fprintf(1, 'preparing mask selection... ')
    %COORDINATES ARE ROTATED 21degrees about the x-axis
    alpha = deg_alpha*pi/180;
    beta = 0*pi/180;
    gamma = 0*pi/180;
    Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
    Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
    Rz = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];

    Rx = Rx*Ry*Rz;
    im_cortex_name = [subjectName,'cortex.jpg'];

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
    rot_edgeL = (Rx*total_edge')';
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
    
    %rotate left edge
    rot_left_edgeL = (Rx*left_edge')';
    
    %rotate right edge
    rot_right_edgeL = (Rx*right_edge')';
    
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
    
    %ROTATE CORTEX LOCATIONS TO MATCH TOP DOWN LOOK
%     rot_cortexL_full = (Rx*cortexL_full')';
    
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
   
    
    %ROTATE everything else too
    rot_censul = (Rx*censul')';


    %IMAGE PROCESSING OF CORTEX PICTURE
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

    %+leftshift, -rightshift
    imA_resize = imA((min(im_i)-3):(max(im_i)+1),(min(im_j)-12):(max(im_j)+16),:);
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
    figure(1)
    set(gcf,'Color',[1 1 1]);
    imshow(imA_resize2,'Xdata',[mincortX maxcortX],'YData',[maxcortY mincortY])
    axis on
    axis([mincortX maxcortX mincortY maxcortY])
    axis xy
    hold on;
    %scatter3(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),7,'y');
    %scatter3(loc_cortex(:,1),loc_cortex(:,2),zeros(length(loc_cortex(:,1)),1),7,'m')
    plot3(rot_censul(:,1),rot_censul(:,2),rot_censul(:,3),'.','Markersize',15,'color','b')
    %scatter3(rot_censul(:,1),rot_censul(:,2),rot_censul(:,3),3,'g');
    %scatter3(rot_targetL(:,1),rot_targetL(:,2),rot_targetL(:,3),7,'r');

    %plot midline
    %line([midline1_x midline1_x],[mincortY maxcortY],'Color','red','Linestyle','--')
    %plot(midline1_x,y_locs,'color','m','Linewidth',2,'Linestyle','--')
    %plot(midline2_x,y_locs,'color','m','Linewidth',2,'Linestyle','--')

    view(0,90)
    axis on
    axis equal
    axis([mincortX maxcortX mincortY maxcortY])

    %resize imA_resize2 so that pixel size = samp_inc
    imA_resize3 = imresize(imA_resize2,[length(bigy),length(bigx)],'bilinear');
    imA_resize4 = imresize(imA_resize2,[length(bigy)-1,length(bigx)-1],'bilinear');

    %done with preparing mask selection
    fprintf(1, 'done\n')
    
    %CREATE MASK OF ROTATED TARGET LOCATIONS
%     if (createmaskflag)
    fprintf(1, 'create mask... ')

    %create mask array- limits of ROTATED target locations                    
    mask = zeros(length(rot_y),length(rot_x));
    [mask_m,mask_n] = size(mask);
    edge_mask = mask;
    right_edge_mask = edge_mask;
    left_edge_mask = edge_mask;
    censul_mask = edge_mask;

    Xmat = repmat(rot_x,mask_m,1);
    Ymat = repmat(rot_y,mask_n,1)';

    %fill mask with locations of targets
    for loc_ind=1:length(rot_targetL(:,1))
        %dx = rot_targetL(loc_ind,1); 
        %dy = rot_targetL(loc_ind,2);
        [minval,minrow] = min((rot_targetL(loc_ind,1)-Xmat).^2 + (rot_targetL(loc_ind,2)-Ymat).^2);
        [minval,mincol] = min(minval);

        loc_y = minrow(mincol); 
        loc_x = mincol;
        mask(loc_y,loc_x) = 1;
    end

    %fill edge_mask with locations of total edges
    for loc_ind=1:length(rot_edgeL(:,1))
        [minval,minrow] = min((rot_edgeL(loc_ind,1)-Xmat).^2 + (rot_edgeL(loc_ind,2)-Ymat).^2);
        [minval,mincol] = min(minval);

        loc_y = minrow(mincol); 
        loc_x = mincol;
        edge_mask(loc_y,loc_x) = 1;
    end
    
    %fill right_edge_mask with locations of right edge
    for loc_ind=1:length(rot_right_edgeL(:,1))
        [minval,minrow] = min((rot_right_edgeL(loc_ind,1)-Xmat).^2 + (rot_right_edgeL(loc_ind,2)-Ymat).^2);
        [minval,mincol] = min(minval);

        loc_y = minrow(mincol); 
        loc_x = mincol;
        right_edge_mask(loc_y,loc_x) = 1;
    end
    
    %fill left_edge_mask with locations of left edge
    for loc_ind=1:length(rot_left_edgeL(:,1))
        [minval,minrow] = min((rot_left_edgeL(loc_ind,1)-Xmat).^2 + (rot_left_edgeL(loc_ind,2)-Ymat).^2);
        [minval,mincol] = min(minval);

        loc_y = minrow(mincol); 
        loc_x = mincol;
        left_edge_mask(loc_y,loc_x) = 2;
        mask(loc_y,loc_x) = 2;
    end
    
    %fill censul_mask with locations of central sulcus
    for loc_ind=1:length(rot_censul(:,1))
        [minval,minrow] = min((rot_censul(loc_ind,1)-Xmat).^2 + (rot_censul(loc_ind,2)-Ymat).^2);
        [minval,mincol] = min(minval);

        loc_y = minrow(mincol); 
        loc_x = mincol;
        censul_mask(loc_y,loc_x) = 1;
    end

    %fill masks
    mask2 = imclose(mask,se);
    mask3 = imfill(mask2,'holes');
    mask4 = imopen(mask3,se);
    mask5 = mask4+2*censul_mask;   %show proposed pts with original edge pts selected in Curry


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
    bigcensul_mask = bigmask5;

    yind_lt_corner_mask5 = round(l_bigy*(maxcortY-maxrot_Y)/(maxcortY-mincortY));
    xind_lt_corner_mask5 = max(round(l_bigx*(minrot_X-mincortX)/(maxcortX-mincortX)),1);
    
    bigmask5(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(mask5);
    bigmask4(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(mask4);
    bigmask(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(mask);

    bigedge_mask(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(edge_mask);
    bigcensul_mask(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(censul_mask);

    %load old bigmask
    if ~isempty(allmasks)
        bigmask4 = flipud(allmasks);
    end
    
    hold on;
    h2 = imagesc(bigmask5,[0 3]);
    colormap hsv
    set(h2,'AlphaData',.25);
    set(gca,'CLim',[0 9]);
    axis tight

%     xlimits = get(gca,'XLim');
%     ylimits = get(gca,'YLim');
%     h3 = imagesc(imA_resize3);
    
    %----------------------------------------------------------
    %User picks mask
    but = 1;
    mask_value = 1;
    while (but >= 1 & but <=3) | (but >= 48 & but <=60 ) | (but >= 186 & but <=187 )
        [xi,yi,but] = ginput(1);
        if (but >= 1 & but <=3) | (but >= 48 & but <=60)
            if but == 49    %button 1 on keyboard
                mask_value = 1; %1 signifies right side, M1
            elseif but == 50    %button 2 on keyboard
                mask_value = 2; %2 signifies left side, M1
            elseif but == 51    %button 3 on keyboard
                mask_value = 3; %3 signifies right side, S1
            elseif but == 52    %button 4 on keyboard
                mask_value = 4; %4 signifies left side, S1
            elseif but == 53    %button 5 on keyboard
                mask_value = 5; %5 signifies right side, PM
            elseif but == 54    %button 6 on keyboard
                mask_value = 6; %6 signifies left side, PM
            elseif but == 55    %button 7 on keyboard
                mask_value = 7; %7 signifies right side, SMA
            elseif but == 56    %button 8 on keyboard
                mask_value = 8; %8 signifies left side, SMA
            elseif but == 57    %button 5 on keyboard
                mask_value = 9; %9 signifies right side, preSMA
            elseif but == 48    %button 6 on keyboard
                mask_value = 10; %10 signifies left side, preSMA
            elseif but == 186    %button 7 on keyboard
                mask_value = 11; %11 signifies left side, Pa
            elseif but == 187    %button 8 on keyboard
                mask_value = 12; %12 signifies left side, Pa
                
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
                set(h2,'CData',bigmask4+3*bigedge_mask+0*bigcensul_mask)
            end
        end
    end
    allmasks = flipud(bigmask4);
    %close figure
    close(100)

    %done with mask creation
    fprintf(1, 'done\n')       
    %---------------------------------------------------
    
    bigmask4 = allmasks;
    mask_perim = bwperim(bigmask4,8);
%     mask_perim_nums = union(mask_perim_nums,[areanum]);
    mask_perims = mask_perim;

    %-----------------------------------

    %12e. Make new edge and find new targetL
    %locations

    fprintf(1, 'create new mask... ')

    %CONVERT MASKS to LOCATION MASKS
    %all regions
    [mask_m,mask_n] = find(bigmask4>0);
    for mask_ind = 1:length(mask_m)
        rot_new_edge(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    

    %------------------------------------
    %find individual regions locations 
    for i_area = 1: length (AreaNames)
        %using Curry 5.0 Carolina's pts
        cur_area = AreaNames {i_area};
        
    end
    %find right M1
    [mask_m,mask_n] = find(bigmask4 == 1);
    for mask_ind = 1:length(mask_m)
        right_M1_ROI(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    
    %find left M1
    %even numbers above 0
    [mask_m,mask_n] = find(bigmask4 == 2);
    for mask_ind = 1:length(mask_m)
        left_M1_ROI(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    
    %find right S1
    %even numbers above 0
    [mask_m,mask_n] = find(bigmask4 == 3);
    right_S1_ROI = [];
    for mask_ind = 1:length(mask_m)
        right_S1_ROI(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    
    %find left S1
    %even numbers above 0
    [mask_m,mask_n] = find(bigmask4 == 4);
    left_S1_ROI = [];
    for mask_ind = 1:length(mask_m)
        left_S1_ROI(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    
    %find right PM
    %even numbers above 0
    [mask_m,mask_n] = find(bigmask4 == 5);
    right_PM_ROI = [];
    for mask_ind = 1:length(mask_m)
        right_PM_ROI(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    
    %find left PM
    %even numbers above 0
    [mask_m,mask_n] = find(bigmask4 == 6);
    left_PM_ROI = [];
    for mask_ind = 1:length(mask_m)
        left_PM_ROI(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    
    %find right SMA
    %even numbers above 0
    [mask_m,mask_n] = find(bigmask4 == 7);
    right_SMA_ROI = [];
    for mask_ind = 1:length(mask_m)
        right_SMA_ROI(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    
    %find left SMA
    %even numbers above 0
    [mask_m,mask_n] = find(bigmask4 == 8);
    left_SMA_ROI = [];
    for mask_ind = 1:length(mask_m)
        left_SMA_ROI(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    
     %find right preSMA
    %even numbers above 0
    [mask_m,mask_n] = find(bigmask4 == 9);
    right_preSMA_ROI = [];
    for mask_ind = 1:length(mask_m)
        right_preSMA_ROI(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    
    %find left preSMA
    %even numbers above 0
    [mask_m,mask_n] = find(bigmask4 == 10);
    left_preSMA_ROI = [];
    for mask_ind = 1:length(mask_m)
        left_preSMA_ROI(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
   
      %find right Pa
    %even numbers above 0
    [mask_m,mask_n] = find(bigmask4 == 11);
    right_Pa_ROI = [];
    for mask_ind = 1:length(mask_m)
        right_Pa_ROI(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    
    %find left Pa
    %even numbers above 0
    [mask_m,mask_n] = find(bigmask4 == 12);
    left_Pa_ROI = [];
    for mask_ind = 1:length(mask_m)
        left_Pa_ROI(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
  
    
    %------------------------------------
    
    %find mask perimeter
    [mask_m,mask_n] = find(mask_perim>0);
    for mask_ind = 1:length(mask_m)
        rot_new_perim(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end

    
    regions = rot_new_perim;
    alledges = rot_new_edge;
    fprintf(1, 'done\n')
    
    ROI=rot_new_edge;
%     left_ROI = rot_new_left_edge;
%     right_ROI = rot_new_right_edge;
%     top_cortex_ind=indzs;
    
    
    %----------------------------------------------------------------------
    %go back to curry locations and delete the extra points picked up using matlab--
    %by jun on 09/27/06
    %--------Select the points inside the edge---------------
    fprintf(1, 'delete the extra points... ')
    ROI_ind=[];
    [n,tmp]=size(ROI);

    for i=1:n
        dis2total_edge = sqrt(sum((rot_edgeL-repmat(ROI(i,:),size(rot_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<100 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<25
                ROI_ind = [ROI_ind;ind];
            end
            %------------------------
        end
    end
    ROI_ind = unique(ROI_ind);
    
    
    for i_area = 1: length(AreaNames)
    end
%     %left side
%     left_ROI_ind=[];
%     [left_n,tmp]=size(left_ROI);
% 
%     for i=1:left_n
%         dis2total_edge = sqrt(sum((rot_left_edgeL-repmat(left_ROI(i,:),size(rot_left_edgeL,1),1)).^2,2));
%         if min(dis2total_edge)<100 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
%             %-----go back to cortex index
%             dis2cortex= sqrt(sum((rot_cortexL-repmat(left_ROI(i,:),size(rot_cortexL,1),1)).^2,2));
% 
%             [yind,ind] = min(dis2cortex);
%             if dis2cortex(ind)<25
%                 if i==1
%                     left_ROI_ind = ind;
%                 end
%                 left_ROI_ind = [left_ROI_ind;ind];
%             end
%             %------------------------
%         end
%     end
%     left_ROI_ind = unique(left_ROI_ind);
%     
%     %right side
%     right_ROI_ind=[];
%     [right_n,tmp]=size(right_ROI);
% 
%     for i=1:right_n
%         dis2total_edge = sqrt(sum((rot_right_edgeL-repmat(right_ROI(i,:),size(rot_right_edgeL,1),1)).^2,2));
%         if min(dis2total_edge)<7 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
%             %-----go back to cortex index
%             dis2cortex= sqrt(sum((rot_cortexL-repmat(right_ROI(i,:),size(rot_cortexL,1),1)).^2,2));
% 
%             [yind,ind] = min(dis2cortex);
%             if dis2cortex(ind)<5
%                 if i==1
%                     right_ROI_ind = ind;
%                 end
%                 right_ROI_ind = [right_ROI_ind;ind];
%             end
%             %------------------------
%         end
%     end
%     right_ROI_ind = unique(right_ROI_ind);
    
    %-------------------------------------------------------
    %right M1
    right_M1_ROI_ind=[];
    [right_n,tmp]=size(right_M1_ROI);

    for i=1:right_n
        dis2total_edge = sqrt(sum((rot_right_edgeL-repmat(right_M1_ROI(i,:),size(rot_right_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<7 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(right_M1_ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<5
                if i==1
                    right_M1_ROI_ind = ind;
                end
                right_M1_ROI_ind = [right_M1_ROI_ind;ind];
            end
            %------------------------
        end
    end
    right_M1_ROI_ind = unique(right_M1_ROI_ind);
    
    %left M1
    left_M1_ROI_ind=[];
    [left_n,tmp]=size(left_M1_ROI);

    for i=1:left_n
        dis2total_edge = sqrt(sum((rot_left_edgeL-repmat(left_M1_ROI(i,:),size(rot_left_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<7 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(left_M1_ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<5
                if i==1
                    left_M1_ROI_ind = ind;
                end
                left_M1_ROI_ind = [left_M1_ROI_ind;ind];
            end
            %------------------------
        end
    end
    left_M1_ROI_ind = unique(left_M1_ROI_ind);
    
    %-------------------
    %right S1
    right_S1_ROI_ind=[];
    [right_n,tmp]=size(right_S1_ROI);

    for i=1:right_n
        dis2total_edge = sqrt(sum((rot_right_edgeL-repmat(right_S1_ROI(i,:),size(rot_right_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<7 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(right_S1_ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<5
                if i==1
                    right_S1_ROI_ind = ind;
                end
                right_S1_ROI_ind = [right_S1_ROI_ind;ind];
            end
            %------------------------
        end
    end
    right_S1_ROI_ind = unique(right_S1_ROI_ind);
    
    %left S1
    left_S1_ROI_ind=[];
    [left_n,tmp]=size(left_S1_ROI);

    for i=1:left_n
        dis2total_edge = sqrt(sum((rot_left_edgeL-repmat(left_S1_ROI(i,:),size(rot_left_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<7 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(left_S1_ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<5
                if i==1
                    left_S1_ROI_ind = ind;
                end
                left_S1_ROI_ind = [left_S1_ROI_ind;ind];
            end
            %------------------------
        end
    end
    left_S1_ROI_ind = unique(left_S1_ROI_ind);
    
    %-------------------
    %right PM
    right_PM_ROI_ind=[];
    [right_n,tmp]=size(right_PM_ROI);

    for i=1:right_n
        dis2total_edge = sqrt(sum((rot_right_edgeL-repmat(right_PM_ROI(i,:),size(rot_right_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<7 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(right_PM_ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<5
                if i==1
                    right_PM_ROI_ind = ind;
                end
                right_PM_ROI_ind = [right_PM_ROI_ind;ind];
            end
            %------------------------
        end
    end
    right_PM_ROI_ind = unique(right_PM_ROI_ind);
    
    %left PM
    left_PM_ROI_ind=[];
    [left_n,tmp]=size(left_PM_ROI);

    for i=1:left_n
        dis2total_edge = sqrt(sum((rot_left_edgeL-repmat(left_PM_ROI(i,:),size(rot_left_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<7 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(left_PM_ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<5
                if i==1
                    left_PM_ROI_ind = ind;
                end
                left_PM_ROI_ind = [left_PM_ROI_ind;ind];
            end
            %------------------------
        end
    end
    left_PM_ROI_ind = unique(left_PM_ROI_ind);
    
    %-------------------
    %right SMA
    right_SMA_ROI_ind=[];
    [right_n,tmp]=size(right_SMA_ROI);

    for i=1:right_n
        dis2total_edge = sqrt(sum((rot_right_edgeL-repmat(right_SMA_ROI(i,:),size(rot_right_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<7 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(right_SMA_ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<5
                if i==1
                    right_SMA_ROI_ind = ind;
                end
                right_SMA_ROI_ind = [right_SMA_ROI_ind;ind];
            end
            %------------------------
        end
    end
    right_SMA_ROI_ind = unique(right_SMA_ROI_ind);
    
    %left SMA
    left_SMA_ROI_ind=[];
    [left_n,tmp]=size(left_SMA_ROI);

    for i=1:left_n
        dis2total_edge = sqrt(sum((rot_left_edgeL-repmat(left_SMA_ROI(i,:),size(rot_left_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<7 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(left_SMA_ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<5
                if i==1
                    left_SMA_ROI_ind = ind;
                end
                left_SMA_ROI_ind = [left_SMA_ROI_ind;ind];
            end
            %------------------------
        end
    end
    left_SMA_ROI_ind = unique(left_SMA_ROI_ind);
    %-------------------------------------------------------
    
    %get unique elements for each region
    intersect_ROI_ind = intersect(left_ROI_ind,right_ROI_ind);
    left_ROI_ind = setdiff(left_ROI_ind,intersect_ROI_ind);
    right_ROI_ind = setdiff(right_ROI_ind,intersect_ROI_ind);
    
    %get unique elements for subregions
    left_M1_ROI_ind = intersect(left_ROI_ind,left_M1_ROI_ind);
    left_S1_ROI_ind = intersect(left_ROI_ind,left_S1_ROI_ind);
    left_PM_ROI_ind = intersect(left_ROI_ind,left_PM_ROI_ind);
    left_SMA_ROI_ind = intersect(left_ROI_ind,left_SMA_ROI_ind);
    right_M1_ROI_ind = intersect(right_ROI_ind,right_M1_ROI_ind);
    right_S1_ROI_ind = intersect(right_ROI_ind,right_S1_ROI_ind);
    right_PM_ROI_ind = intersect(right_ROI_ind,right_PM_ROI_ind);
    right_SMA_ROI_ind = intersect(right_ROI_ind,right_SMA_ROI_ind);
    
    %--------------------------------------
    %left side
    
    %give intersection of M1 and S1 to M1
    intersect_ROI_ind = intersect(left_S1_ROI_ind,left_M1_ROI_ind);
    left_S1_ROI_ind = setdiff(left_S1_ROI_ind,intersect_ROI_ind);
    
    %give intersection of S1 and SMA to S1
    intersect_ROI_ind = intersect(left_S1_ROI_ind,left_SMA_ROI_ind);
    left_SMA_ROI_ind = setdiff(left_SMA_ROI_ind,intersect_ROI_ind);
    
    %give intersection of S1 and PM to S1
    intersect_ROI_ind = intersect(left_S1_ROI_ind,left_PM_ROI_ind);
    left_PM_ROI_ind = setdiff(left_PM_ROI_ind,intersect_ROI_ind);
    
    %give intersection of SMA and PM to SMA
    intersect_ROI_ind = intersect(left_SMA_ROI_ind,left_PM_ROI_ind);
    left_PM_ROI_ind = setdiff(left_PM_ROI_ind,intersect_ROI_ind);
    
    %give intersection of SMA and M1 to SMA
    intersect_ROI_ind = intersect(left_SMA_ROI_ind,left_M1_ROI_ind);
    left_M1_ROI_ind = setdiff(left_M1_ROI_ind,intersect_ROI_ind);
    
    %give intersection of PM and M1 to PM
    intersect_ROI_ind = intersect(left_PM_ROI_ind,left_M1_ROI_ind);
    left_M1_ROI_ind = setdiff(left_M1_ROI_ind,intersect_ROI_ind);
    
    %--------------------------------------
    %right side
    
    %give intersection of M1 and S1 to M1
    intersect_ROI_ind = intersect(right_S1_ROI_ind,right_M1_ROI_ind);
    right_S1_ROI_ind = setdiff(right_S1_ROI_ind,intersect_ROI_ind);
    
    %give intersection of S1 and SMA to S1
    intersect_ROI_ind = intersect(right_S1_ROI_ind,right_SMA_ROI_ind);
    right_SMA_ROI_ind = setdiff(right_SMA_ROI_ind,intersect_ROI_ind);
    
    %give intersection of S1 and PM to S1
    intersect_ROI_ind = intersect(right_S1_ROI_ind,right_PM_ROI_ind);
    right_PM_ROI_ind = setdiff(right_PM_ROI_ind,intersect_ROI_ind);
    
    %give intersection of SMA and PM to SMA
    intersect_ROI_ind = intersect(right_SMA_ROI_ind,right_PM_ROI_ind);
    right_PM_ROI_ind = setdiff(right_PM_ROI_ind,intersect_ROI_ind);
    
    %give intersection of SMA and M1 to SMA
    intersect_ROI_ind = intersect(right_SMA_ROI_ind,right_M1_ROI_ind);
    right_M1_ROI_ind = setdiff(right_M1_ROI_ind,intersect_ROI_ind);
    
    %give intersection of PM and M1 to PM
    intersect_ROI_ind = intersect(right_PM_ROI_ind,right_M1_ROI_ind);
    right_M1_ROI_ind = setdiff(right_M1_ROI_ind,intersect_ROI_ind);
    
    %------------------------------------------------
    
    %zeros refer to right ROI, 1's refer to left ROI
    %ROI_ind = [left_ROI_ind;right_ROI_ind];
    ROI_ind2 = [left_M1_ROI_ind;left_S1_ROI_ind;left_PM_ROI_ind;left_SMA_ROI_ind;right_M1_ROI_ind;right_S1_ROI_ind;right_PM_ROI_ind;right_SMA_ROI_ind];
    ROI_ind = ROI_ind2;
    
    ROI_ind_group = [ones(length(left_ROI_ind),1);zeros(length(right_ROI_ind),1)];
    ROI_ind_group2 = [(ones(length(left_M1_ROI_ind),1)+1);...
        (ones(length(left_S1_ROI_ind),1)+3);...
        (ones(length(left_PM_ROI_ind),1)+5);...
        (ones(length(left_SMA_ROI_ind),1)+7);...
        (zeros(length(right_M1_ROI_ind),1)+1);...
        (zeros(length(right_S1_ROI_ind),1)+3);...
        (zeros(length(right_PM_ROI_ind),1)+5);...
        (zeros(length(right_SMA_ROI_ind),1)+7)];
    
    if length(ROI_ind_group) == length(ROI_ind_group2)
        fprintf(1,'ROI_ind_group indices match ROI_ind_group2\n');
    end
    
    if yj_debug==1
        ROI_cortex=rot_cortexL(ROI_ind,:);
        left_ROI_cortex=rot_cortexL(left_ROI_ind,:);
        right_ROI_cortex=rot_cortexL(right_ROI_ind,:);
        
        left_M1_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==2),:);
        left_S1_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==4),:);
        left_PM_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==6),:);
        left_SMA_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==8),:);
        right_M1_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==1),:);
        right_S1_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==3),:);
        right_PM_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==5),:);
        right_SMA_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==7),:);
        
        %3D plot of points
        figure(2)
%         plot3(rot_cortexL_full(:,1),rot_cortexL_full(:,2),rot_cortexL_full(:,3),'.', 'MarkerSize',3,'Color','g');
%         hold on;
        plot3(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),'.', 'MarkerSize',3);
        hold on
        plot3(ROI(:,1),ROI(:,2),ROI(:,3),'m.', 'MarkerSize',3);
        plot3(right_ROI_cortex(:,1),right_ROI_cortex(:,2),right_ROI_cortex(:,3),'r.', 'MarkerSize',6);
        plot3(left_ROI_cortex(:,1),left_ROI_cortex(:,2),left_ROI_cortex(:,3),'g.', 'MarkerSize',6);
        
        %2D plot of points on top of cortex image
        figure(1)
        plot3(right_ROI_cortex(:,1),right_ROI_cortex(:,2),right_ROI_cortex(:,3),'r.', 'MarkerSize',6);
        plot3(left_ROI_cortex(:,1),left_ROI_cortex(:,2),left_ROI_cortex(:,3),'g.', 'MarkerSize',6);
        
        plot3(right_M1_ROI_cortex(:,1),right_M1_ROI_cortex(:,2),right_M1_ROI_cortex(:,3),'r.', 'MarkerSize',8);
        plot3(right_S1_ROI_cortex(:,1),right_S1_ROI_cortex(:,2),right_S1_ROI_cortex(:,3),'m.', 'MarkerSize',8);
        plot3(right_PM_ROI_cortex(:,1),right_PM_ROI_cortex(:,2),right_PM_ROI_cortex(:,3),'c.', 'MarkerSize',8);
        plot3(right_SMA_ROI_cortex(:,1),right_SMA_ROI_cortex(:,2),right_SMA_ROI_cortex(:,3),'b.', 'MarkerSize',8);
        plot3(left_M1_ROI_cortex(:,1),left_M1_ROI_cortex(:,2),left_M1_ROI_cortex(:,3),'g.', 'MarkerSize',8);        
        plot3(left_S1_ROI_cortex(:,1),left_S1_ROI_cortex(:,2),left_S1_ROI_cortex(:,3),'y.', 'MarkerSize',8);        
        plot3(left_PM_ROI_cortex(:,1),left_PM_ROI_cortex(:,2),left_PM_ROI_cortex(:,3),'c.', 'MarkerSize',8);        
        plot3(left_SMA_ROI_cortex(:,1),left_SMA_ROI_cortex(:,2),left_SMA_ROI_cortex(:,3),'k.', 'MarkerSize',8);        
        
        
        
    end
    fprintf(1, 'done\n')
    fprintf(1, 'Please check the results by rotating the figure...\n')
    reply = input('Do you want to try removing some points? Y/N [Y]: ','s');
    while isempty(reply)
        reply = input('Do you want to try removing some points? Y/N [Y]: ','s');
    end

    while (reply == 'Y' | reply=='y')
        
        reply2 = [];
        while isempty(reply2)
            reply2 = input('Automatic, Manual, or Quit? A/M/Q [A]: ','s');
        end
        
        if ~isequal(reply2,'A') && ~isequal(reply2,'a') && ~isequal(reply2,'M') && ~isequal(reply2,'m') && ~isequal(reply2,'Q') && ~isequal(reply2,'q')
            reply2 = 'A';
        end
        
        clear ind
        
        if reply2 == 'A' || reply2 == 'a'
            
            s_maxclust_num = [];
            while isempty(s_maxclust_num)
                s_maxclust = input('How many clusters are there? [try max of 4]: ','s');
                if isempty(s_maxclust)
                    s_maxclust = '4';
                end
                s_maxclust_num = str2num(s_maxclust);
                if isempty(s_maxclust_num)
                    printf(1,'Wrong input\n');
                end
            end

            fprintf(1, 'Displaying clusters... ')
            ROI=rot_cortexL(ROI_ind,:);
            Dis=pdist(ROI);
            Z = linkage(Dis);
            T = cluster(Z,'maxclust',str2num(s_maxclust));
            for clust_no=1:str2num(s_maxclust)
                ind{clust_no}=find(T==clust_no);
                clust_length(clust_no)=length(ind{clust_no});
            end
            [max_clust,biggest_clust]=max(clust_length);
            [min_clust,smallest_clust] = min(clust_length);
            %         ROI_ind=ROI_ind(ind{biggest_clust});
            
            figure(4)
            clf
            hold on;
            for clust_no=1:str2num(s_maxclust)
                ROI_cluster=rot_cortexL(ROI_ind(ind{clust_no}),:);

                if clust_no ~= biggest_clust
                    color = rand([3 1]);
                    marker_size = 12;
                else
                    color = 'r';
                    marker_size = 5;
                end

                if length(ROI_cluster)>=3
                    plot3(ROI_cluster(:,1),ROI_cluster(:,2),ROI_cluster(:,3),'.','Color',color,'MarkerSize',marker_size);
                end
            end
            
            reply = input('Do you want to remove these points? Y/N [Y]: ','s');
            if isempty(reply) | reply == 'Y' | reply=='y'
                ROI_ind=ROI_ind(ind{biggest_clust});
                ROI_ind_group = ROI_ind_group(ind{biggest_clust});
                ROI_ind_group2 = ROI_ind_group2(ind{biggest_clust});
            end
            
        elseif reply2 == 'M' || reply2 == 'm'
            %use ginput to manually select point to delete
            figure(4)
            clf
            hold on;
            
            ROI_cluster=rot_cortexL(ROI_ind,:);
            color = 'r';
            marker_size = 5;
            
            plot3(ROI_cluster(:,1),ROI_cluster(:,2),ROI_cluster(:,3),'.','Color',color,'MarkerSize',marker_size);
            
            plot3(ROI_cluster(ROI_ind_group>0,1),ROI_cluster(ROI_ind_group>0,2),ROI_cluster(ROI_ind_group>0,3),'.','Color','g','MarkerSize',marker_size+1);
            
            a1 = gca;
            %a2 = axes;
            %set(a2,'Position',get(a1,'Position'));
            plot3(ROI_cluster(ROI_ind_group2==1,1),ROI_cluster(ROI_ind_group2==1,2),ROI_cluster(ROI_ind_group2==1,3),'.','Color','r','MarkerSize',marker_size+1);
            plot3(ROI_cluster(ROI_ind_group2==2,1),ROI_cluster(ROI_ind_group2==2,2),ROI_cluster(ROI_ind_group2==2,3),'.','Color','m','MarkerSize',marker_size+1);
            
            if ~isempty(ROI_cluster(ROI_ind_group2==3,1))
                plot3(ROI_cluster(ROI_ind_group2==3,1),ROI_cluster(ROI_ind_group2==3,2),ROI_cluster(ROI_ind_group2==3,3),'.','Color','g','MarkerSize',marker_size+1);
            end
            if ~isempty(ROI_cluster(ROI_ind_group2==4,1))
                plot3(ROI_cluster(ROI_ind_group2==4,1),ROI_cluster(ROI_ind_group2==4,2),ROI_cluster(ROI_ind_group2==4,3),'.','Color','y','MarkerSize',marker_size+1);
            end
            if ~isempty(ROI_cluster(ROI_ind_group2==5,1))
                plot3(ROI_cluster(ROI_ind_group2==5,1),ROI_cluster(ROI_ind_group2==5,2),ROI_cluster(ROI_ind_group2==5,3),'.','Color','c','MarkerSize',marker_size+1);
            end
            if ~isempty(ROI_cluster(ROI_ind_group2==6,1))
                plot3(ROI_cluster(ROI_ind_group2==6,1),ROI_cluster(ROI_ind_group2==6,2),ROI_cluster(ROI_ind_group2==6,3),'.','Color','b','MarkerSize',marker_size+1);
            end
            if ~isempty(ROI_cluster(ROI_ind_group2==7,1))
                plot3(ROI_cluster(ROI_ind_group2==7,1),ROI_cluster(ROI_ind_group2==7,2),ROI_cluster(ROI_ind_group2==7,3),'.','Color','b','MarkerSize',marker_size+1);
            end
            if ~isempty(ROI_cluster(ROI_ind_group2==8,1))
                plot3(ROI_cluster(ROI_ind_group2==8,1),ROI_cluster(ROI_ind_group2==8,2),ROI_cluster(ROI_ind_group2==8,3),'.','Color','c','MarkerSize',marker_size+1);
            end
            %axes(a1);
            
            biggest_clust = 1;
            ind{1} = 1:length(ROI_cluster);
            
            delete_pt = 1;
            del_pts = [];
            
            % del1 = plot3(ROI_cluster(:,1),ROI_cluster(:,2),ROI_cluster(:,3),'.','Color',color,'MarkerSize',marker_size,'Visible','off');
            del1 = plot3(ROI_cluster(ROI_ind_group2>0,1),ROI_cluster(ROI_ind_group2>0,2),ROI_cluster(ROI_ind_group2>0,3),'.','Color',color,'MarkerSize',marker_size,'Visible','off');
            
            %view x-z view
            v = [1 0 0 -0.5; 0 0 1 -0.5; 0 1 0 8.1603; 0 0 0 1];
            view(v);
            % del2 = plot3(ROI_cluster(:,1),ROI_cluster(:,2),ROI_cluster(:,3),'.','Color',color,'MarkerSize',marker_size);
            but = 1;
            while but == 1 | but == 3 | but == 2 | but == 49 | but == 50
                fprintf(1,'paused for rotation (press enter to go on)...\n')
                pause
                [xi,yi,but] = ginput(1);
                [p v i_min_pt] = select3d;
                %set(del2,'XData',ROI_cluster(vi,1),'YData',ROI_cluster(vi,2),'ZData',ROI_cluster(vi,3),'Linestyle','.','Color','g','MarkerSize',12);
                
                if but==1 | but == 3 | but == 2 | but == 49 | but == 50
                    if but == 49
                        delete_pt = 1;  %change to delete mode
                    elseif but == 50
                        delete_pt = 0;  %change to undelete mode
                    else
                        %find pt closest to mouse click and select it- use select3d.m file
                        if ~isempty(v)
                            [y_min_pt,i_min_pt] = min((ROI_cluster(:,1)-v(1)).^2 + (ROI_cluster(:,2)-v(2)).^2 + (ROI_cluster(:,3)-v(3)).^2);
                            i_min_pt = i_min_pt(1);
                            % ROI_cluster(i_min_pt,:);

                            if delete_pt
                                ind{1} = setdiff(ind{1},i_min_pt);
                                del_pts = sort(unique([del_pts i_min_pt]));

                                set(del1,'XData',ROI_cluster(del_pts,1),'YData',ROI_cluster(del_pts,2),'ZData',ROI_cluster(del_pts,3),'Linestyle','.','Color','k','MarkerSize',12,'Visible','on');
                            else
                                ind{1} = sort(union(ind{1},i_min_pt));
                                del_pts = sort(setdiff(del_pts,i_min_pt));

                                set(del1,'XData',ROI_cluster(del_pts,1),'YData',ROI_cluster(del_pts,2),'ZData',ROI_cluster(del_pts,3),'Linestyle','.','Color','k','MarkerSize',12,'Visible','on');
                            end
                        end
                    end
                end
            end
            
            reply3 = [];
            while isempty(reply3) 
                reply3 = input('Do you want to remove these points? Y/N [Y]: ','s');
            end
            if reply3 == 'Y' | reply3 =='y'
                ROI_ind=ROI_ind(ind{biggest_clust});
                ROI_ind_group = ROI_ind_group(ind{biggest_clust});
                ROI_ind_group2 = ROI_ind_group2(ind{biggest_clust});
            end
        else
            reply = 'q';
        end
            
        
        
        if yj_debug==1
            ROI_cortex=rot_cortexL(ROI_ind,:);
            left_ROI_cortex = rot_cortexL(ROI_ind(ROI_ind_group>0),:);
            right_ROI_cortex = rot_cortexL(ROI_ind(ROI_ind_group<1),:);
            
            figure(2)
            clf
            plot3(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),'.', 'MarkerSize',3);
            hold on
%             plot3(left_ROI_cortex(:,1),left_ROI_cortex(:,2),left_ROI_cortex(:,3),'g.', 'MarkerSize',6);
%             plot3(right_ROI_cortex(:,1),right_ROI_cortex(:,2),right_ROI_cortex(:,3),'r.', 'MarkerSize',6);
            plot3(ROI_cluster(ROI_ind_group2==1,1),ROI_cluster(ROI_ind_group2==1,2),ROI_cluster(ROI_ind_group2==1,3),'.','Color','r','MarkerSize',marker_size+1);
            plot3(ROI_cluster(ROI_ind_group2==2,1),ROI_cluster(ROI_ind_group2==2,2),ROI_cluster(ROI_ind_group2==2,3),'.','Color','m','MarkerSize',marker_size+1);
            
            if ~isempty(ROI_cluster(ROI_ind_group2==3,1))
                plot3(ROI_cluster(ROI_ind_group2==3,1),ROI_cluster(ROI_ind_group2==3,2),ROI_cluster(ROI_ind_group2==3,3),'.','Color','g','MarkerSize',marker_size+1);
            end
            if ~isempty(ROI_cluster(ROI_ind_group2==4,1))
                plot3(ROI_cluster(ROI_ind_group2==4,1),ROI_cluster(ROI_ind_group2==4,2),ROI_cluster(ROI_ind_group2==4,3),'.','Color','y','MarkerSize',marker_size+1);
            end
            if ~isempty(ROI_cluster(ROI_ind_group2==5,1))
                plot3(ROI_cluster(ROI_ind_group2==5,1),ROI_cluster(ROI_ind_group2==5,2),ROI_cluster(ROI_ind_group2==5,3),'.','Color','c','MarkerSize',marker_size+1);
            end
            if ~isempty(ROI_cluster(ROI_ind_group2==6,1))
                plot3(ROI_cluster(ROI_ind_group2==6,1),ROI_cluster(ROI_ind_group2==6,2),ROI_cluster(ROI_ind_group2==6,3),'.','Color','b','MarkerSize',marker_size+1);
            end
            if ~isempty(ROI_cluster(ROI_ind_group2==7,1))
                plot3(ROI_cluster(ROI_ind_group2==7,1),ROI_cluster(ROI_ind_group2==7,2),ROI_cluster(ROI_ind_group2==7,3),'.','Color','b','MarkerSize',marker_size+1);
            end
            if ~isempty(ROI_cluster(ROI_ind_group2==8,1))
                plot3(ROI_cluster(ROI_ind_group2==8,1),ROI_cluster(ROI_ind_group2==8,2),ROI_cluster(ROI_ind_group2==8,3),'.','Color','c','MarkerSize',marker_size+1);
            end
        end

    end
    
    if yj_debug==1
        ROI_cortex=rot_cortexL(ROI_ind,:);
        left_ROI_cortex = rot_cortexL(ROI_ind(ROI_ind_group>0),:);
        right_ROI_cortex = rot_cortexL(ROI_ind(ROI_ind_group<1),:);
        
        left_M1_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==2),:);
        left_S1_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==4),:);
        left_PM_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==6),:);
        left_SMA_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==8),:);
        right_M1_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==1),:);
        right_S1_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==3),:);
        right_PM_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==5),:);
        right_SMA_ROI_cortex = rot_cortexL(ROI_ind2(ROI_ind_group2==7),:);
        
        figure(5)
        set(gcf,'Color',[1 1 1]);
        imshow(imA_resize2,'XData',[mincortX maxcortX],'YData',[maxcortY mincortY])
        axis on
        axis([mincortX maxcortX mincortY maxcortY])
        axis xy
        hold on;

        plot3(ROI_cluster(ROI_ind_group2==1,1),ROI_cluster(ROI_ind_group2==1,2),ROI_cluster(ROI_ind_group2==1,3),'.','Color','r','MarkerSize',marker_size+1);
        plot3(ROI_cluster(ROI_ind_group2==2,1),ROI_cluster(ROI_ind_group2==2,2),ROI_cluster(ROI_ind_group2==2,3),'.','Color','m','MarkerSize',marker_size+1);

        if ~isempty(ROI_cluster(ROI_ind_group2==3,1))
            plot3(ROI_cluster(ROI_ind_group2==3,1),ROI_cluster(ROI_ind_group2==3,2),ROI_cluster(ROI_ind_group2==3,3),'.','Color','g','MarkerSize',marker_size+1);
        end
        if ~isempty(ROI_cluster(ROI_ind_group2==4,1))
            plot3(ROI_cluster(ROI_ind_group2==4,1),ROI_cluster(ROI_ind_group2==4,2),ROI_cluster(ROI_ind_group2==4,3),'.','Color','y','MarkerSize',marker_size+1);
        end
        if ~isempty(ROI_cluster(ROI_ind_group2==5,1))
            plot3(ROI_cluster(ROI_ind_group2==5,1),ROI_cluster(ROI_ind_group2==5,2),ROI_cluster(ROI_ind_group2==5,3),'.','Color','c','MarkerSize',marker_size+1);
        end
        if ~isempty(ROI_cluster(ROI_ind_group2==6,1))
            plot3(ROI_cluster(ROI_ind_group2==6,1),ROI_cluster(ROI_ind_group2==6,2),ROI_cluster(ROI_ind_group2==6,3),'.','Color','b','MarkerSize',marker_size+1);
        end
        if ~isempty(ROI_cluster(ROI_ind_group2==7,1))
            plot3(ROI_cluster(ROI_ind_group2==7,1),ROI_cluster(ROI_ind_group2==7,2),ROI_cluster(ROI_ind_group2==7,3),'.','Color','b','MarkerSize',marker_size+1);
        end
        if ~isempty(ROI_cluster(ROI_ind_group2==8,1))
            plot3(ROI_cluster(ROI_ind_group2==8,1),ROI_cluster(ROI_ind_group2==8,2),ROI_cluster(ROI_ind_group2==8,3),'.','Color','c','MarkerSize',marker_size+1);
        end
        view(0,90)
        axis on
        axis equal
        axis([mincortX maxcortX mincortY maxcortY])
        title([subjectName])
    end
    save ([Base_dir,subjectName,'\Points\','Mask.mat'],'ROI_ind','ROI_ind_group','ROI_ind_group2','Rx','allmasks');
end