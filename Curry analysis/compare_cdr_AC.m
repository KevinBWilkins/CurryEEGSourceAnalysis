
function compare_cdr_AC(varargin)

%control_mat_list = {'GMrobot_RE_new','YZrobot1_RE_new','BSrobot1_RE_new','JDrobot_RE'};
control_mat_list = {'GMrobot_RE_new'};
%control_mat_list = {'GMrobot_0RE','YZrobot1_0RE_new','BSrobot1_0RE','JDrobot_0RE'};
%control_mat_list = {'GMrobot_25RE','YZrobot2_25RE','BSrobot1_25RE','JDrobot_25RE'};

stroke_mat_list = {'NTrobot1_RE_new'};
stroke_mat_list = {};
%stroke_mat_list = {'NTrobot1_0RE'};
%stroke_mat_list = {'NTrobot2_25RE_new'};
stroke_hand_list = [1]; %right = 1, left=0;

%------------------------
numControl = length(control_mat_list);
numStroke = length(stroke_mat_list);
numSubjects = numControl+numStroke;

for i = 1:numSubjects
    if i<=numControl
        m1 = control_mat_list{i};
    else
        m1 = stroke_mat_list{i-numControl};
    end
    
    matname = ['C:\Albert Chen\matlab\Curry analysis\cdr\', m1,'_cdr.mat'];
    load(matname)
    %variables in mat file = cdr, Subject
    
    subjectName = cdr.subjectName;
    plotstart = cdr.plotstart;
    plotstop = cdr.plotstop;
    newtargetL = cdr.newtargetL;
    norm_newtargetL = cdr.norm_newtargetL;
    newtargetC = cdr.newtargetC;
    maxtargetC = cdr.maxtargetC;
    newThreshold = cdr.newThreshold;
    norm_mean_locs = cdr.norm_mean_locs;
    norm_loc_cortex2 = cdr.norm_loc_cortex2;
    imA_resize3 = cdr.imA_resize3;
    rot_new_perim = cdr.rot_new_perim;
    norm_new_perim = cdr.norm_new_perim;
    rot_censul = cdr.rot_censul;
    norm_rot_censul = cdr.norm_rot_censul;
    norm_newrot_XI = cdr.norm_newrot_XI;
    norm_newrot_YI = cdr.norm_newrot_YI;
    newrot_ZI = cdr.newrot_ZI;
    newmaxrot_Z = cdr.newmaxrot_Z;
    v_caxis = cdr.v_caxis;
    new_mask_perim_filled_resize = cdr.new_mask_perim_filled_resize;
    
    %pick plot color
    if strfind(matname,'25RE')
        m_color = 'm';
    elseif strfind(matname,'0RE')
        m_color = 'g';
    else
        m_color = 'b';
    end
    
    %------------------------------------
    %reproduce figure from Curry_analysis_AC
    figure(i)
    hold on;
    %top-down image of cortex
    imshow(imA_resize3,'XData',[-1 1],'YData',[1 -1])
    axis on
    axis([-1 1 -1 1])
    hold on;

    plot3(norm_new_perim(:,1),norm_new_perim(:,2),rot_new_perim(:,3),'LineWidth',7,'color','r','Marker','.','Linestyle','.')
    scatter3(norm_loc_cortex2(:,1),norm_loc_cortex2(:,2),zeros(length(norm_loc_cortex2(:,1)),1),7,'m')
    plot3(norm_rot_censul(:,1),norm_rot_censul(:,2),rot_censul(:,3),'Markersize',15,'color','b','LineStyle','.')
    plot3(norm_newtargetL(:,1),norm_newtargetL(:,2),newtargetL(:,3),'LineWidth',7,'color','y','LineStyle','.')
    plot3(norm_mean_locs(1,:),norm_mean_locs(2,:),newmaxrot_Z*ones(1,length(norm_mean_locs(1,:))),'LineWidth',7,'color',m_color,'LineStyle','.')

    CI = mean(newtargetC(:,plotstart:plotstop),2);
    norm_CI = griddata(norm_newtargetL(:,1),norm_newtargetL(:,2),CI,norm_newrot_XI,norm_newrot_YI);
    h3 = surface('XData',norm_newrot_XI,'YData',norm_newrot_YI,'ZData',newrot_ZI,'CData',norm_CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6);
    set(h3,'FaceAlpha','flat','AlphaData',new_mask_perim_filled_resize);
    set(gca,'ALim',[0 1]);

    view(0,90)
    axis equal
    axis xy
    caxis(v_caxis);
    %axis on;axis equal
    axis([-1 1 -1 1])
    title(m1)
    close(i)
    
    %put all brains on one figure, with transparency equal to 1/number of people
    
    %align y coordinates by mean of y locations
    mean_y = mean(norm_newtargetL(:,2));
    norm_newtargetL(:,2) = norm_newtargetL(:,2)-mean_y;
    norm_new_perim(:,2) = norm_new_perim(:,2)-mean_y;
    norm_rot_censul(:,2) = norm_rot_censul(:,2)-mean_y;
    norm_mean_locs(2,:) = norm_mean_locs(2,:)-mean_y;
    
    if i<=numControl
        figure(10);
    else
        figure(20);
    end
    
    hold on;
    %top-down image of cortex
    i1 = imshow(imA_resize3,'XData',[-1 1],'YData',[1 -1]);
    hold on;
    set(i1,'AlphaData',1/numSubjects);
    axis on
    view(0,90)
    axis equal
    axis xy
    axis([-1 1 -1 1])
    
    if i<=numControl
        figure(11);
    else
        figure(21);
    end
    hold on;
    p1 = plot3(norm_new_perim(:,1),norm_new_perim(:,2),zeros(length(rot_new_perim(:,3)),1),'LineWidth',7,'color','r','Marker','.','Linestyle','.');
    %set(p1,'AlphaData',1/numSubjects);
    hold on;
    plot3(norm_rot_censul(:,1),norm_rot_censul(:,2),zeros(length(rot_censul(:,3)),1),'Markersize',15,'color','b','LineStyle','.')
    hold on;
    plot3(norm_newtargetL(:,1),norm_newtargetL(:,2),zeros(length(newtargetL(:,3)),1),'LineWidth',7,'color','y','LineStyle','.')
    hold on;
    plot3(norm_mean_locs(1,:),norm_mean_locs(2,:),newmaxrot_Z*ones(1,length(norm_mean_locs(1,:))),'LineWidth',7,'color',m_color,'LineStyle','.')
    hold on;
    axis on
    view(0,90)
    axis equal
    axis xy
    axis([-1 1 -1 1])
    
    %calculate average COG
    
    
    
    
    %perform normalization transformation
    
    %transformation matrix
    
    
    
    %------------------------------------------------
end