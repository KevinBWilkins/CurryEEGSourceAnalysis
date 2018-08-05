% function [censul,targetL,targetC_BL,targetC,XI,YI,ZI]=select_ROI (cortexL,targetL,edge)                
function [ROI_ind,ROI_ind_group,Rx,allmasks]=select_ROI_AC_S1(varargin) 
% This function returns the rotated cortex locations and the index of the region of interests. The input
% varargin can be:
%    1. subjectName
%    2. deg_alpha
%    3. remakemask
%    4. SubjectGroup

% subjectName = 'WSrobot';
% deg_alpha = -21;
% exmaple:
% [ROI_ind,ROI_ind_group,Rx,allmasks]=select_ROI_AC_S1('DW',-21,[],'control');

warning off MATLAB:griddata:DuplicateDataPoints

subjectName = varargin{1};
deg_alpha = varargin{2};
remakemask = 1;
SubjectGroup = varargin{4};

timePoints=1;
areanum=1;
yj_debug=1;

right_edge = [];
left_edge = [];
allmasks = [];
Base_dir=['F:\data\inverse_results\',SubjectGroup,'\'];

if exist ([Base_dir,subjectName,'\Points\','Mask.mat'],'file')
%     createmaskflag = 1;
    allmasks = [];
    ROI_ind_group = [];
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
    %using Curry 5.0 Carolina's pts
    cdr_file=['Ne_Lelbow3'];
%     cdr_BL_file=['RE_new_mri_BL'];
    
    disp('********************');
    disp(cdr_file)
%     cdr_file_name=[Base_dir,subjectName,'\AKpre_CC\',cdr_file,'.cdr'];
    cdr_file_name=[Base_dir,subjectName,'\new\',cdr_file,'.cdr'];
    
%     subjectNames_allregion = ['FNrobot1', 'FNrobot2', 'YZrobot1', 'LS1', 'LS2','LSpost1','LSpost2','GMrobot','JS2','CM103006','JS2new','AKpost1','AKpost2','AKpre1','AKpre2'];
    subjectNames_newregion = ['WSrobot', 'WSrobot2', 'WSrobot3','JWrobot','JWrobot2','JDrobot','BSrobot1'];
    subjectNames_amputee = ['AK'];

    subjectNames_allregion = [];
%     subjectNames_newregion = [];
    %EDGE IS TAKEN HERE TO MEAN "REGION"
    if strfind(subjectNames_allregion,subjectName)
        edge1_file_name=[Base_dir,subjectName,'\Points\','allregions_ac_new.pom'];%all regions
        [edge1,Ecount,ENR]=read_Curry_file4_AC(edge1_file_name,'LOCATION',0,0);
        total_edge = [edge1];
        right_edge = total_edge;
        left_edge = total_edge;
        censul_file_name1=[Base_dir,subjectName,'\Points\','CS_cc.pom'];   %central sulcus
%     elseif strfind(subjectNames_allregion,subjectName)
%         edge1_file_name=[Base_dir,subjectName,'\Points\','S1_rt_cc.pom'];%S1_rt
%         edge2_file_name=[Base_dir,subjectName,'\Points\','S1_lt_cc.pom'];%S1_lt
% 
%                 %READ IN LOCATIONS OF POINTS
%         [edge1,Ecount,ENR]=read_Curry_file4_AC(edge1_file_name,'LOCATION',0,0);
%         [edge2,Ecount,ENR]=read_Curry_file4_AC(edge2_file_name,'LOCATION',0,0);
% 
% %         total_edge = [edge1;edge2;edge3;edge4];
%         total_edge = [edge1;edge2];
% 
%         right_edge = [edge1];
%         left_edge = [edge2];
% 
%         censul_file_name1=[Base_dir,subjectName,'\Points\','CS_cc.pom'];   %central sulcus
    else
        edge1_file_name=[Base_dir,subjectName,'\Points\','S1_rt_cc.pom'];%S1_rt
        edge2_file_name=[Base_dir,subjectName,'\Points\','S1_lt_cc.pom'];%S1_lt
        
%         if strfind(subjectNames_newregion,subjectName)
%             edge9_file_name=[Base_dir,subjectName,'\Points\','new_region.pom'];%newregion
% %             edge9_file_name=[Base_dir,subjectName,'\Points\','SMA_new.pom'];%newregion
%             %[edge9,Ecount,ENR]=read_Curry_file4_AC(edge9_file_name,'LOCATION',0,0);
%             [edge9,Ecount,LNR,TM9]=read_Curry_file3_TM(edge9_file_name,'LOCATION',0,0);
%         end
        
        %READ IN LOCATIONS OF POINTS
        [edge1,Ecount,ENR]=read_Curry_file4_AC(edge1_file_name,'LOCATION',0,0);
        [edge2,Ecount,ENR]=read_Curry_file4_AC(edge2_file_name,'LOCATION',0,0);
        
        
        total_edge = [edge1;edge2];
        right_edge = [edge1];
        left_edge = [edge2];
%         if strfind(subjectNames_newregion,subjectName)
%             
%             %multiply by transformation matrix
%             TM2 = TM9(2:4,1);    %translation
%             TM1 = TM9(2:4,2:4);  %rotation and scaling
%                     
%             newedge9 = (TM1*edge9')' + repmat(TM2',size(edge9,1),1);
%             
%             %make sure newregion has same transformation matrix as SMA_rt_cc
%             TM2 = TM7(2:4,1);    %translation
%             TM1 = TM7(2:4,2:4);  %rotation and scaling
%                     
%             newedge9 = (TM1*edge9')' + repmat(TM2',size(edge9,1),1);
%             
%             
%             total_edge = [edge1;edge2;edge3;edge4;edge5;edge6;edge7;edge8;newedge9];
%             right_edge = [edge1;edge3;edge5;edge7;newedge9];
%             left_edge = [edge2;edge4;edge6;edge8;newedge9];
%             
%         
%         end
        %Read in central sulcus file, may be different name
        censul_file_name1=[Base_dir,subjectName,'\Points\','CS_cc.pom'];   %central sulcus
    end


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
       
    %if do not want to apply transformation
    cortexL = cortexL_mri_voxel;
    
    %--------------------------------------------
    %get all cortex locations corresponding closest to edge points
    fprintf(1, 'get ROI locations... ')
    %use matlab picked ROI regions for total_edge if you want to just modify an old mask
    if exist ([Base_dir,subjectName,'\Points\','Mask.mat'],'file')
        
        %comment out for updating new ROI
%         total_edge = cortexL(ROI_ind,:);
%         targetL = total_edge;
%         
%         if ~isempty(ROI_ind_group)
%             left_ROI_ind = ROI_ind(ROI_ind_group>0);
%             left_targetL = cortexL(left_ROI_ind,:);
%             right_ROI_ind = ROI_ind(ROI_ind_group<1);
%             right_targetL = cortexL(right_ROI_ind,:);
%             
%             %1's belong to left side, 0's belong to right side
%             targetL_group = ismember(total_edge,left_targetL);
%             left_group = targetL_group>0;
%             left_group = left_group(:,1);
%             left_group = left_group>0;
%             right_group = ~left_group;
%             
%             %update left and right edges
%             left_edge = targetL(left_group,:);
%             right_edge = targetL(right_group,:);
%         end
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

    %get highest z-values at each x-y coordinate (only top surface of cortex)//??? Why?-Jun
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
    imshow([mincortX maxcortX],[maxcortY mincortY],imA_resize2)
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
    %imshow(imA_resize2)
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
    colormap vga
    set(h2,'AlphaData',.25);

    but = 1;
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
                set(h2,'CData',bigmask4+0*bigedge_mask+0*bigcensul_mask)
            end
        end
    end
    allmasks = flipud(bigmask4);
    %close figure
    close(100)

    %done with mask creation
    fprintf(1, 'done\n')       

    bigmask4 = allmasks;
    mask_perim = bwperim(bigmask4,8);
%     mask_perim_nums = union(mask_perim_nums,[areanum]);
    mask_perims = mask_perim;

    %-----------------------------------

    %12e. Make new edge and find new targetL
    %locations

    fprintf(1, 'create new mask... ')

    %CONVERT MASKS to LOCATION MASKS
    [mask_m,mask_n] = find(bigmask4>0);
    for mask_ind = 1:length(mask_m)
        rot_new_edge(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    
    %find locations of left side
    [mask_m,mask_n] = find(bigmask4>1);
    for mask_ind = 1:length(mask_m)
        rot_new_left_edge(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end
    
    %find locations of right side
    [mask_m,mask_n] = find(bigmask4<2 & bigmask4>0);
    for mask_ind = 1:length(mask_m)
        rot_new_right_edge(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end

    [mask_m,mask_n] = find(mask_perim>0);
    for mask_ind = 1:length(mask_m)
        rot_new_perim(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end

    regions = rot_new_perim;
    alledges = rot_new_edge;
    fprintf(1, 'done\n')
    
    ROI=rot_new_edge;
    left_ROI = rot_new_left_edge;
    right_ROI = rot_new_right_edge;
%     top_cortex_ind=indzs;
    
    
    %-------go back to curry locations and delete the extra points picked up using matlab--by jun on
    %09/27/06--------------------------------------------------------
    %--------Select the points inside the edge---------------
    fprintf(1, 'delete the extra points... ')
    ROI_ind=[];
    [n,tmp]=size(ROI);

    for i=1:n
        dis2total_edge =  sqrt(sum((rot_edgeL-repmat(ROI(i,:),size(rot_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<200 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<5
                ROI_ind = [ROI_ind;ind];
            end
            %------------------------
        end
    end
    ROI_ind = unique(ROI_ind);
    
    %left side
    left_ROI_ind=[];
    [left_n,tmp]=size(left_ROI);

    for i=1:left_n
        dis2total_edge =  sqrt(sum((rot_left_edgeL-repmat(left_ROI(i,:),size(rot_left_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<7 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(left_ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<5
                if i==1
                    left_ROI_ind = ind;
                end
                left_ROI_ind = [left_ROI_ind;ind];
            end
            %------------------------
        end
    end
    left_ROI_ind = unique(left_ROI_ind);
    
    %right side
    right_ROI_ind=[];
    [right_n,tmp]=size(right_ROI);

    for i=1:right_n
        dis2total_edge =  sqrt(sum((rot_right_edgeL-repmat(right_ROI(i,:),size(rot_right_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<7 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(right_ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<5
                if i==1
                    right_ROI_ind = ind;
                end
                right_ROI_ind = [right_ROI_ind;ind];
            end
            %------------------------
        end
    end
    right_ROI_ind = unique(right_ROI_ind);
    
    intersect_ROI_ind = intersect(left_ROI_ind,right_ROI_ind);
    left_ROI_ind = setdiff(left_ROI_ind,intersect_ROI_ind);
    right_ROI_ind = setdiff(right_ROI_ind,intersect_ROI_ind);
    
    %zeros refer to right ROI, 1's refer to left ROI
    ROI_ind = [left_ROI_ind;right_ROI_ind];
    %ROI_ind_group = ismember(ROI_ind,left_ROI_ind);
    ROI_ind_group = [ones(length(left_ROI_ind),1);zeros(length(right_ROI_ind),1)];
    
    if yj_debug==1
        ROI_cortex=rot_cortexL(ROI_ind,:);
        left_ROI_cortex=rot_cortexL(left_ROI_ind,:);
        right_ROI_cortex=rot_cortexL(right_ROI_ind,:);
        figure(2)
%         plot3(rot_cortexL_full(:,1),rot_cortexL_full(:,2),rot_cortexL_full(:,3),'.', 'MarkerSize',3,'Color','g');
%         hold on;
        plot3(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),'.', 'MarkerSize',3);
        hold on
        plot3(ROI(:,1),ROI(:,2),ROI(:,3),'m.', 'MarkerSize',3);
        plot3(right_ROI_cortex(:,1),right_ROI_cortex(:,2),right_ROI_cortex(:,3),'r.', 'MarkerSize',6);
        plot3(left_ROI_cortex(:,1),left_ROI_cortex(:,2),left_ROI_cortex(:,3),'g.', 'MarkerSize',6);
        
        figure(1)
        plot3(right_ROI_cortex(:,1),right_ROI_cortex(:,2),right_ROI_cortex(:,3),'r.', 'MarkerSize',6);
        plot3(left_ROI_cortex(:,1),left_ROI_cortex(:,2),left_ROI_cortex(:,3),'g.', 'MarkerSize',6);
        
        
    end
    fprintf(1, 'done\n')
    fprintf(1, 'Please check the results by rotating the figure...\n')
    reply = input('Do you want to try removing some points? Y/N [Y]: ','s');
    if isempty(reply)
        reply = 'Y';
    end

    while (isempty(reply) | reply == 'Y' | reply=='y')
        
        reply2 = input('Automatic, Manual, or Quit? A/M/Q [A]: ','s');
        if ~isequal(reply2,'A') && ~isequal(reply2,'a') && ~isequal(reply2,'M') && ~isequal(reply2,'m') && ~isequal(reply2,'Q') && ~isequal(reply2,'q')
            reply2 = 'A';
        end
        
        clear ind
        
        if reply2 == 'A' || reply2 == 'a'
            
            s_maxclust_num = [];
            clust_length = [];
            clear len_small2big_clust small2big_clust
            while isempty(s_maxclust_num)
                s_maxclust = input('How many clusters are there? [try max of 15]: ','s');
                if isempty(s_maxclust)
                    s_maxclust = '15';
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
            
                   
 %--------Jun's new update---------------------                   
%             for clust_no=1:str2num(s_maxclust)
%                ind{clust_no}=find(T==clust_no); 
%                clust_length(clust_no)=length(ind{clust_no});
%             end
            [len_small2big_clust,small2big_clust]=sort(clust_length);
            len_small2big_clust
            clust_included=1;
            for clust_included_no=1:clust_no
                if (sum(len_small2big_clust([end-clust_included_no+1:end]))<0.8*sum(len_small2big_clust))
                    clust_included=clust_included+1;
                else
                    break;
                end
            end
    %         if (len_small2big_clust(end)>0.9*sum(len_small2big_clust))
    %             ROI_ind=ROI_ind(ind{small2big_clust(end)});
    %         else 
                %fprintf(1, 'too many clusters? ... ');
            tmp=ind([small2big_clust([end-clust_included+1:end])]);

            tmp1=cell2mat(tmp');
%             ROI_ind=ROI_ind(tmp1);
%             s_maxclust = num2str(str2num(s_maxclust)-1);

%------------display the remaining cluster using red, and other clusters using other colors-------------------------            
			figure(4)
			clf
			hold on;
			for clust_no=1:str2num(s_maxclust)
                current_clust_ind = small2big_clust (end - clust_no +1);
                current_clust_cortexL =  rot_cortexL(ROI_ind(ind{current_clust_ind}),:);
			%                 ROI_cluster=rot_cortexL(ROI_ind(ind{clust_no}),:);
			
                if clust_no > clust_included
                    color = rand([3 1]);
                    marker_size = 12;
                else
                    color = 'r';
                    marker_size = 5;
                end
			
                if length(current_clust_cortexL)>=3
                    plot3(current_clust_cortexL(:,1),current_clust_cortexL(:,2),current_clust_cortexL(:,3),'.','Color',color,'MarkerSize',marker_size);
                end
			end
            reply = input('Do you want to remove these points? Y/N [Y]: ','s');
            if isempty(reply) | reply == 'Y' | reply=='y'
                ROI_ind=ROI_ind(tmp1);
                ROI_ind_group = ROI_ind_group(tmp1);
            end
            
        elseif reply2 == 'M' || reply2 == 'm'
            %use ginput to manually select point to delete
            figure(4)
            clf
            hold on;
            
            ROI_cluster=rot_cortexL(ROI_ind,:);
            color = 'r';
            marker_size = 5;
            % plot3(ROI_cluster(ROI_ind_group<1,1),ROI_cluster(ROI_ind_group<1,2),ROI_cluster(ROI_ind_group<1,3),'.','Color','r','MarkerSize',marker_size);

            plot3(ROI_cluster(:,1),ROI_cluster(:,2),ROI_cluster(:,3),'.','Color',color,'MarkerSize',marker_size);
            plot3(ROI_cluster(ROI_ind_group>0,1),ROI_cluster(ROI_ind_group>0,2),ROI_cluster(ROI_ind_group>0,3),'.','Color','g','MarkerSize',marker_size+1);
            biggest_clust = 1;
            ind{1} = 1:length(ROI_cluster);
            
            delete_pt = 1;
            del_pts = [];
            del1 = plot3(ROI_cluster(:,1),ROI_cluster(:,2),ROI_cluster(:,3),'.','Color',color,'MarkerSize',marker_size,'Visible','off');

            %view x-z view
            v = [1 0 0 -0.5; 0 0 1 -0.5; 0 1 0 8.1603; 0 0 0 1];
%             view(v);
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
%                         [y_min_pt,i_min_pt] = min((ROI_cluster(:,1)-xi).^2 + (ROI_cluster(:,3)-zi).^2);
%                         i_min_pt = i_min_pt(1);
%                         ROI_cluster(i_min_pt,:)
                        
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
            
            reply = input('Do you want to remove these points? Y/N [Y]: ','s');
            if isempty(reply) | reply == 'Y' | reply=='y'
                ROI_ind=ROI_ind(ind{biggest_clust});
                ROI_ind_group = ROI_ind_group(ind{biggest_clust});
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
            plot3(left_ROI_cortex(:,1),left_ROI_cortex(:,2),left_ROI_cortex(:,3),'g.', 'MarkerSize',6);
            plot3(right_ROI_cortex(:,1),right_ROI_cortex(:,2),right_ROI_cortex(:,3),'r.', 'MarkerSize',6);
            update_disp_points_in_cortex_pic(rot_cortexL, ROI_ind, ROI_ind_group, mincortX, maxcortX, maxcortY, mincortY, imA_resize2);
        end

    end
    
%     if yj_debug==1
% %         ROI_cortex=rot_cortexL(ROI_ind,:);
% %         left_ROI_cortex = rot_cortexL(ROI_ind(ROI_ind_group>0),:);
% %         right_ROI_cortex = rot_cortexL(ROI_ind(ROI_ind_group<1),:);
% %         
% %         figure(5)
% %         set(gcf,'Color',[1 1 1]);
% %         imshow([mincortX maxcortX],[maxcortY mincortY],imA_resize2)
% %         axis on
% %         axis([mincortX maxcortX mincortY maxcortY])
% %         axis xy
% %         hold on;
% % 
% %         plot3(left_ROI_cortex(:,1),left_ROI_cortex(:,2),left_ROI_cortex(:,3),'g.', 'MarkerSize',6);
% %         plot3(right_ROI_cortex(:,1),right_ROI_cortex(:,2),right_ROI_cortex(:,3),'r.', 'MarkerSize',6);
% % 
% %         view(0,90)
% %         axis on
% %         axis equal
% %         axis([mincortX maxcortX mincortY maxcortY])
% %         title([subjectName])
%         update_disp_points_in_cortex_pic(rot_cortexL, ROI_ind, ROI_ind_group, mincortX, maxcortX, maxcortY, mincortY, imA_resize2);
%     end
    save ([Base_dir,subjectName,'\Points\','Mask.mat'],'ROI_ind','ROI_ind_group','Rx','allmasks');
    figure (5)
    figName = [Base_dir,subjectName,'\Points\','Mask-pic.jpg'];
    saveas (gcf, figName , 'jpg');
end