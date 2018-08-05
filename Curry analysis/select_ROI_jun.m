% function [censul,targetL,targetC_BL,targetC,XI,YI,ZI]=select_ROI (cortexL,targetL,edge)                
function [ROI_ind,Rx]=select_ROI (subjectName,deg_alpha) 
% This function returns the rotated cortex locations and the index of the region of interests. The input
% varargin can be two types:
%    1. (cortexL, targetL)
%    2. (subjectName)

timePoints=1;
areanum=1;
yj_debug=1;
Base_dir='F:\data\inverse_results\stroke-post\';

if exist ([Base_dir,subjectName,'\points\','Mask.mat'],'file')
    createmaskflag = 1;
    load ([Base_dir,subjectName,'\points\','Mask.mat'])
else
    createmaskflag=0;
    %----------------------------------------------------------------------
    %using Curry 5.0 Carolina's pts
    cdr_file=['abd.cdr'];
    cdr_BL_file=['abd.cdr'];

    disp('********************');
    disp(cdr_file)
    cdr_file_name=[Base_dir,subjectName,'\new\',cdr_file];
    cdr_BL_file_name=[Base_dir,subjectName,'\new\',cdr_BL_file];
    
    edge1_file_name=[Base_dir,subjectName,'\points\','M1_rt_cc.pom'];%M1_rt
    edge2_file_name=[Base_dir,subjectName,'\points\','M1_lt_cc.pom'];%M1_lt
    edge3_file_name=[Base_dir,subjectName,'\points\','S1_rt_cc.pom'];%S1_rt
    edge4_file_name=[Base_dir,subjectName,'\points\','S1_lt_cc.pom'];%S1_lt
    edge5_file_name=[Base_dir,subjectName,'\points\','PM_rt_cc.pom'];%PM_rt
    edge6_file_name=[Base_dir,subjectName,'\points\','PM_lt_cc.pom'];%PM_lt
    edge7_file_name=[Base_dir,subjectName,'\points\','SMA_rt_cc.pom'];%SMA_rt
    edge8_file_name=[Base_dir,subjectName,'\points\','SMA_lt_cc.pom'];%SMA_lt

    %READ IN LOCATIONS OF POINTS
    [edge1,Ecount,ENR]=read_Curry_file4_AC(edge1_file_name,'LOCATION',0,0);
    [edge2,Ecount,ENR]=read_Curry_file4_AC(edge2_file_name,'LOCATION',0,0);
    [edge3,Ecount,ENR]=read_Curry_file4_AC(edge3_file_name,'LOCATION',0,0);
    [edge4,Ecount,ENR]=read_Curry_file4_AC(edge4_file_name,'LOCATION',0,0);
    [edge5,Ecount,ENR]=read_Curry_file4_AC(edge5_file_name,'LOCATION',0,0);
    [edge6,Ecount,ENR]=read_Curry_file4_AC(edge6_file_name,'LOCATION',0,0);
    [edge7,Ecount,ENR]=read_Curry_file4_AC(edge7_file_name,'LOCATION',0,0);
    [edge8,Ecount,ENR]=read_Curry_file4_AC(edge8_file_name,'LOCATION',0,0);

    total_edge = [edge1;edge2;edge3;edge4;edge5;edge6;edge7;edge8];
    areanum = 1;
    %--------------------------------------------------------------

    %Read in central sulcus file, may be different name
    censul_file_name1=[Base_dir,subjectName,'\points\','CS_cc.pom'];   %central sulcus

    %if only one central sulcus file
    [censul1,Ecount,ENR]=read_Curry_file3(censul_file_name1,'LOCATION',0,0);
    censul2 = censul1;
    censul=[censul1];
    %----------------------------------------------------------------------    
    %get all cortex locations
    fprintf(1, 'get cortex locations... ')
    [cortexL,Lcount,LNR]=read_Curry_file3(cdr_BL_file_name,'LOCATION',0,0);
    fprintf(1, 'done\n')

    %get all cortex strengths corresponding closest to edge points
    fprintf(1, 'get ROI locations using curry points cortex strengths... ')
    [cortexC,targetL]=find_cortex_strengths_AC_jun(cdr_BL_file_name,cortexL,total_edge,timePoints,0);
    %[cortexL,targetL,targetC]=find_target_PC3_AC(cdr_file_name,cortexL,edge,timePoints,0);
    fprintf(1, 'done\n')
    
    %--------------------------------------------




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
    im_cortex_file_name=[Base_dir,subjectName,'\points\',im_cortex_name];%image file location
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
    imshow([mincortX maxcortX],[maxcortY mincortY],imA_resize2)
    axis on
    axis([mincortX maxcortX mincortY maxcortY])
    axis xy
    hold on;
    %scatter3(rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),7,'y');
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
%     if (createmaskflag)
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
    bigmask(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(mask);

    bigedge_mask(yind_lt_corner_mask5:yind_lt_corner_mask5+mask_m-1,xind_lt_corner_mask5:xind_lt_corner_mask5+mask_n-1) = flipud(edge_mask);

    hold on;
    h2 = imagesc(bigmask5,[0 3]);
    colormap vga
    set(h2,'AlphaData',.4);

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
    allmasks = flipud(bigmask4);
    %close figure
    close(100)

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

    [mask_m,mask_n] = find(mask_perim>0);
    for mask_ind = 1:length(mask_m)
        rot_new_perim(mask_ind,:) = [bigx(mask_n(mask_ind)) bigy(mask_m(mask_ind)) bigZI(mask_m(mask_ind),mask_n(mask_ind))];
    end

    regions = rot_new_perim;
    alledges = rot_new_edge;
    fprintf(1, 'done\n')
    
    ROI=rot_new_edge;
%     top_cortex_ind=indzs;
    
    
    %-------go back to curry locations and delete the extra points picked up using matlab--by jun on
    %09/27/06--------------------------------------------------------
    %--------Select the points inside the edge---------------
    fprintf(1, 'delete the extra points... ')
    ROI_ind=[];
    [n,tmp]=size(ROI);

    for i=1:n
        dis2total_edge =  sqrt(sum((rot_edgeL-repmat(ROI(i,:),size(rot_edgeL,1),1)).^2,2));
        if min(dis2total_edge)<25 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
            %-----go back to cortex index
            dis2cortex= sqrt(sum((rot_cortexL-repmat(ROI(i,:),size(rot_cortexL,1),1)).^2,2));

            [yind,ind] = min(dis2cortex);
            if dis2cortex(ind)<20
                if i==1
                    ROI_ind = ind;
                end
                ROI_ind = [ROI_ind;ind];
            end
            %------------------------
        end
    end
        
    ROI_ind = unique(ROI_ind);
    
    find(rot_cortexL(:,3)>40);  
    
    if yj_debug==1
        ROI_cortex=rot_cortexL(ROI_ind,:);
        figure
        plot3 (rot_cortexL(find(rot_cortexL(:,3)>40)  ,1),rot_cortexL(find(rot_cortexL(:,3)>40)  ,2),rot_cortexL(find(rot_cortexL(:,3)>40)  ,3),'.', 'MarkerSize',2);
        hold on
        plot3 (ROI(:,1),ROI(:,2),ROI(:,3),'g.', 'MarkerSize',4);
        plot3 (ROI_cortex(:,1),ROI_cortex(:,2),ROI_cortex(:,3),'r.', 'MarkerSize',6);
        %plot3(rot_edgeL(:,1),rot_edgeL(:,2),rot_edgeL(:,3),'k.','MarkerSize',4);
    end
    fprintf(1, 'done\n')
    fprintf(1, 'please check the results by rotating the figure...')
%     reply = input('Do you want further remove some points? Y/N [Y]: ','s');
%     if isempty(reply)
%         reply = 'Y';
%     end
%     %-----------------------------------------------------------------
%     fprintf(1, 'further delete the extra points... ')
% %     ROI=rot_cortexL(ROI_ind,:);
%     bad_ind=[];
%     [n,tmp]=size(ROI_ind);
%     for i=1:n
%         current_ind=ROI_ind(i);
%         other_ind=setdiff(ROI_ind,current_ind);
%         otherPoints=rot_cortexL(other_ind);
%         currentPoint=rot_cortexL(current_ind);
%         dis2otherROI_Points =  sqrt(sum((otherPoints-repmat(currentPoint,size(otherPoints,1),1)).^2,2));
%         dis2otherP=diff(sort(dis2otherROI_Points));
%         
%         if max(dis2otherP)>1 % if miminum distance to the cc points is less than 10, then count it; otherwise, delete it
%             bad_ind=ROI_ind(i);
%         end
%     end
%     ROI_ind = setdiff(ROI_ind,bad_ind);

    %-----------------------------------------------------------------
    reply = input('Do you want further remove some points? Y/N [Y]: ','s');
    if isempty(reply)
        reply = 'Y';
    end

    while (reply == 'Y' | reply=='y')
        clear ind clust_length
        s_maxclust = input('How many clusters there are? [15]: ','s');
        if isempty(s_maxclust)
            s_maxclust = '15';
        end

        fprintf(1, 'further delete the extra points... ')
        ROI=rot_cortexL(ROI_ind,:);
        Dis=pdist(ROI);
        Z = linkage(Dis);
        T = cluster(Z,'maxclust',str2num(s_maxclust));
        for clust_no=1:str2num(s_maxclust)
           ind{clust_no}=find(T==clust_no); 
           clust_length(clust_no)=length(ind{clust_no});
        end
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
        ROI_ind=ROI_ind(tmp1);
        s_maxclust = num2str(str2num(s_maxclust)-1);
%         end
        
        if yj_debug==1
            ROI_cortex=rot_cortexL(ROI_ind,:);
            figure
            plot3 (rot_cortexL(:,1),rot_cortexL(:,2),rot_cortexL(:,3),'.', 'MarkerSize',3);
            hold on
            plot3 (ROI(:,1),ROI(:,2),ROI(:,3),'m.', 'MarkerSize',3);
            plot3 (ROI_cortex(:,1),ROI_cortex(:,2),ROI_cortex(:,3),'r.', 'MarkerSize',6);
        end
        reply = input('Do you want further remove some points? Y/N [Y]: ','s');
        if isempty(reply)
            reply = 'Y';
        end

    end
    
    unRx = [1 0 0; 0 cos(-alpha) sin(-alpha); 0 -sin(-alpha) cos(-alpha)];
    unRy = [cos(-beta) 0 -sin(-beta); 0 1 0; sin(-beta) 0 cos(-beta)];
    unRz = [cos(-gamma) sin(-gamma) 0; -sin(-gamma) cos(-gamma) 0; 0 0 1];
    unRx = unRx*unRy*unRz;
    unrot_newtargetL = (unRx*ROI_cortex')';
    write_new_pom_file(Base_dir,subjectName,unrot_newtargetL);
    
    save ([Base_dir,subjectName,'\points\','Mask.mat'],'ROI_ind','Rx' );
end

% unRx = [1 0 0; 0 cos(-1*alpha) sin(-1*alpha); 0 -sin(-1*alpha) cos(-1*alpha)];
% unRy = [cos(-1*beta) 0 -sin(-1*beta); 0 1 0; sin(-1*beta) 0 cos(-1*beta)];
% unRz = [cos(-1*gamma) sin(-1*gamma) 0; -sin(-1*gamma) cos(-1*gamma) 0; 0 0 1];
% unRx = unRx*unRy*unRz;
% unrot_newtargetL = (unRx*ROI_cortex')';
% write_new_pom_file(Base_dir,subjectName,unrot_newtargetL);


% [censul,targetL,targetC_BL,targetC,XI,YI,ZI]

function write_new_pom_file(Base_dir,subjectName,unrot_newtargetL)

num_locs = length(unrot_newtargetL(:,1));
fid1 = fopen([Base_dir,subjectName,'\points\','new_region','.pom'],'w');

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

