function plot_data_on_cortex(subjectName, deg_alpha, Base_dir, targetL, censul, targetC, cortexL, index, taskColorMap, firstFlag)
% this function plot the inverse calculation result on a cortex.
% Inputs: subjectName: a string showing subject name
%         deg_alpha: rotation angle, usually equals to -21
%         Base_dir: base directory of the cdr file
%         targetL: target location, a matrix with location number x 3
%         censul: locations of central sulcus
%         targetC: current density on target locations, include data in one
%         time window. If want to get a clean target picture without
%         activities lower than threshould, the input should be all the
%         points will values larget than threshould
%         cortexL: all the locations of cortex
%         index: the index of the number of window in time domain


    %1. COORDINATES ARE ROTATED 21degrees about the x-axis
        alpha = deg_alpha*pi/180;
        beta = 0*pi/180;
        gamma = 0*pi/180;
        Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
        Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
        Rz = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];

        Rx = Rx*Ry*Rz;
    
    %2. Rotate cortex and target points
        rot_cortexL=(Rx*cortexL')';
        rot_targetL=(Rx*targetL')';
        rot_censul = (Rx*censul')';
    
    

    %3. Generate CORTEX LOCATIONS in a grid for plot TO MATCH TOP DOWN LOOK
        samp_inc = 0.8;
        mincortX=min(rot_cortexL(:,1));
        maxcortX=max(rot_cortexL(:,1));
        mincortY=min(rot_cortexL(:,2));
        maxcortY=max(rot_cortexL(:,2));
        bigx=mincortX:samp_inc:maxcortX;
        bigy=mincortY:samp_inc:maxcortY;
        [bigXI,bigYI] = meshgrid(bigx,bigy);
        l_bigy = length(bigy);      %dimensions of big cortex locations
        l_bigx = length(bigx);

    %4. find limits of rotated target locations

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
        
    %5. plot the cortex.jpg file first
    
        %IMAGE PROCESSING OF CORTEX PICTURE
        %used to overlay region of interest onto picture of cortex
        if firstFlag == 1
            im_cortex_name = [subjectName,'cortex.jpg'];
            im_cortex_file_name=[Base_dir,subjectName,'\Points\',im_cortex_name];%image file location

        elseif firstFlag == 2
            im_cortex_name = [subjectName,'cortex+abd.jpg'];
            im_cortex_file_name=[Base_dir,subjectName,'\pic\',im_cortex_name];%image file location
       
        end
        fig_h=figure(102);
        if firstFlag < 3
            imA = imread(im_cortex_file_name);       %raw picture
            imB = im2bw(imA,.95);               %make black and white
            se = strel('square',5);             
            imC = imclose(imB,se);              %close up holes
            imD = ~imC;                         %flip white and black
            imE = imfill(imD,'holes');          %fill holes
            [im_i,im_j] = find(imE>0);          %find filled area
            imF = imE(min(im_i):max(im_i),min(im_j):max(im_j));     %
            imG = bwperim(imF,8);

            leftshift = 0;
            rightshift = 0;
            imA_resize = imA((min(im_i)-3):(max(im_i)+1),(min(im_j)-12+leftshift):(max(im_j)+16-rightshift),:);
        %             imA_resize = imA((min(im_i)-3):(max(im_i)+1),(min(im_j)-16):(max(im_j)+17),:);
        %             imA_resize = imA(min(im_i):max(im_i),min(im_j):max(im_j),:);
            %fit cortex picture to cortex pt locations

            [cortex_m,cortex_n] = size(imG);
            imA_resize2 = imresize(imA_resize,[cortex_m,cortex_n],'bilinear');


        %6. normalize the current strength on the target. This is a general
        %normalization. We also can do a norm within each window.
            [locationNo,TimeNo]=size(targetC);
    %         targetC_sum=sum(targetC,2);
    %         targetC_norm=targetC_sum/max(targetC_sum);


        %7. start to plot
             
            set(fig_h,'Color',[1 1 1])
    %         index=0;
    %-------------the following section should be the same main funcion:
    %target_ccd_OI_timebins.m
    %         winOverlap=1;
    %         winSize=13;
    %         timePoints=155;
    %         LoopNo=fix((timePoints-winOverlap)/(winSize-winOverlap));
    %         needed_len=(LoopNo-1)*(winSize-winOverlap)+winSize;
    %         startPoints=timePoints-needed_len-1;


    %         for phTime = 1:LoopNo
    %             index=index+1;
    %             cur_window = (phTime-1)*(winSize-winOverlap)+1+startPoints:(phTime-1)*(winSize-winOverlap)+winSize+startPoints;
                %8 or 16 subplots
        %         subplot(4,4,index)
            subplot(3,4,index)
            imshow(imA_resize2,'XData',[mincortX maxcortX],'YData',[maxcortY mincortY])
            axis([mincortX maxcortX mincortY maxcortY])
            axis off
            hold on

%             timestart_seg = timenow;
%             timestop_seg = timenow + timeinc-1;
%             if timestop_seg > timestop
%                 timestop_seg = timestop;
%             end

            %SURFACE PLOT of active areas, all same amplitude
%             C = mean(targetC(:,:),2);
    %                 [thr_i,thr_j] = find(C<newThreshold);
    %                 Thr_ind = sub2ind(size(C),thr_i,thr_j);
    %                 C_bin = ones(size(C));
    %                 C_bin(Thr_ind) = 0;
            CI_bin = griddata(rot_targetL(:,1),rot_targetL(:,2),targetC,rot_XI,rot_YI);
            surface('XData',rot_XI,'YData',rot_YI,'ZData',rot_ZI,'CData',CI_bin,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.4);
%             contour(rot_XI, rot_YI, CI_bin);
%             set(h3,'FaceAlpha','flat','AlphaData',new_mask_perim_filled_resize);
%             set(gca,'ALim',[0 1]);
            taskColorMap=deblank(taskColorMap);
            colormap (taskColorMap)
        elseif firstFlag == 3
            plot3(rot_targetL(:,1),rot_targetL(:,2),rot_targetL(:,3),'mo');
        end
            view(0,90)
%             TR = timeRange_total{i};
%             title([num2str(TR(1)+(TR(2)-TR(1))*(timestart_seg-1)/(timePoints_total(i)-1),5),' to ',num2str(TR(1)+(TR(2)-TR(1))*(timestop_seg-1)/(timePoints_total(i)-1),5),'ms'],'Fontsize',14,'FontWeight','bold')
            axis off
            caxis([0 1])
            axis equal
            axis xy
            axis tight
%         end
            subject_pic_Dir=[Base_dir,subjectName,'\pic'];
            if ~exist(subject_pic_Dir,'dir')
                mkdir(subject_pic_Dir);
            end
            
            if firstFlag == 1
                fig_file_name=[Base_dir,subjectName,'\pic\',subjectName,'cortex+abd.jpg'];   %index of ROI
            elseif firstFlag == 2
                fig_file_name=[Base_dir,subjectName,'\pic\',subjectName,'cortex+abd+ef.jpg'];   %index of ROI
            else
                fig_file_name=[Base_dir,subjectName,'\pic\',subjectName,'cortex+final.jpg'];   %index of ROI
            end
            saveas (fig_h, fig_file_name, 'jpg');

% 
% %     axis xy
%     hold on;
%     plot3(rot_censul(:,1),rot_censul(:,2),rot_censul(:,3),'.','Markersize',15,'color','b')
%     for phTime=1:LoopNo
%         plot3(CoG(:,phTime,i,1),CoG(:,phTime,i,2),CoG(:,phTime,i,3),'LineStyle','none','Marker','.','MarkerFaceColor','b','EraseMode','none', 'MarkerSize',5,'color',[(LoopNo-phTime)/LoopNo phTime/LoopNo phTime/LoopNo]);
%         drawnow
%         hold on
%     end

    