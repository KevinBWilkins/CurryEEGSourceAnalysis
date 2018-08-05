function [division,division_C] =  divide_loc_current(targetL_all,target_wh_time_ROI,target_area, piont_shift)

step = 3.8*2.7;
yj_color =['y.';'g.';'c.';'y.';'g.';'c.';'y.';'g.';'c.';'y.';'g.';'c.';'y.';'g.';'c.';'y.';'g.';'c.';'y.';'g.';'c.'];

figure (1)
% clf
hold on
for i = 1:length(target_area)
    cur_area = targetL_all{target_area(i)};
    cur_area_C = target_wh_time_ROI {target_area(i)};
    cur_area_C = cur_area_C (: , piont_shift+1:end); %this is only for tim post3.
    if target_area (i) == 2
        cur_area(:,1) = -1*cur_area(:,1);
    end
    [tmp,first_point_ind] = min (cur_area(:,1));
%     [tmp,last_point_ind] = max (cur_area(:,1));
%     total_dis = norm(cur_area(first_point_ind,:)-cur_area(last_point_ind,:));
    
    %------find the middle point in the first line, 
    %definded as the points which is close enough to the middle line--------------------
    begin_ind = first_point_ind;
    cur_line = find_line(begin_ind, cur_area);          
%     mid_point = [mean(cur_line(:,1)),mean(cur_line(:,2)),mean(cur_line(:,3))];
    
    
    %------define the first line function, which passing the first middle point and
    %follow the first line direction---------------------------------------------
    [tmp,highest_point_ind] = max (cur_line(:,2));
    [tmp,lowest_point_ind] = min (cur_line(:,2));
    highest_point = cur_line(highest_point_ind,:);
    lowst_point = cur_line(lowest_point_ind,:);
    clear cur_line tmp
    %-----calculate the distance from all points to the line definded by
    %highest and lowest line
    allP2Line_distance = P2Line_dis(cur_area, highest_point, lowst_point);
%     for i=1:size(cur_area,1)
%         cur_point = cur_area(i,:);
%         allP2Line_distance (i) = norm( cross ((lowst_point-highest_point),(highest_point-cur_point)))/ norm (lowst_point-highest_point);
%     end
    
    cur_division_ind = find (allP2Line_distance<step);
%     division {1} = cur_area (cur_division_ind,:);
%     cur_division = division {1};

    tmp = cur_area (cur_division_ind,:);
    division_C{i,1} = cur_area_C(cur_division_ind,:);
    tmp1=tmp;
    if target_area (i) == 2
        tmp1(:,1) = -1*tmp1(:,1);
    end
    division {i,1} = tmp1;
    
    plot3(tmp1(:,1),tmp1(:,2),tmp1(:,3),'k.','MarkerSize',10);
    clear tmp1
%     pause
    [tmp1,last_point_ind] = max (tmp(:,1));
    cur_line = find_line(last_point_ind, tmp);          
    post_mid_point = [mean(cur_line(:,1)),mean(cur_line(:,2)),mean(cur_line(:,3))];
    cur_area (cur_division_ind,:) =[];
    cur_area_C (cur_division_ind,:) =[]; 
    pre_mid_point = post_mid_point;
    clear tmp1 tmp
    k=1;
    while size(cur_area,1)>10
        k=k+1;
        post_mid_point_tmp = [pre_mid_point(1)+step, pre_mid_point(2), pre_mid_point(3)];
        
%         all_dis = norm(cur_area-repmat(post_mid_point_tmp,size(cur_area,1),1));
        % find the point that is closest to the temp post mid point
        P2Pset_distance = P2Pset_dis (post_mid_point_tmp, cur_area);
        [tmp, post_mid_point_tmp_ind] = min(P2Pset_distance);
        %find the line with x close enought to  the temp post mid point
        cur_line = find_line(post_mid_point_tmp_ind, cur_area);
        post_mid_point_tmp = [mean(cur_line(:,1)),mean(cur_line(:,2)),mean(cur_line(:,3))];
        clear cur_line
        %calculate distance between all the points inside the previous division to the line defined by
        %the two mid points
%         tmp = division {k-1};
        P2L_dis = P2Line_dis(cur_area, pre_mid_point, post_mid_point_tmp);
        % define points on the line
        tmp_ind = P2L_dis<3.8;
        cur_line = cur_area(tmp_ind,:);
        %calculate the distance from all point on the line to the
        %pre_mid_point
%         P2P_dis = norm(cur_line-repmat(pre_mid_point,size(cur_line,1),1));
        P2P_dis = P2Pset_dis (pre_mid_point, cur_line);
        P2P_dis = P2P_dis - step;
        [tmp,tmp_ind] = min(abs(P2P_dis));
        post_mid_point = cur_line (tmp_ind,:);
        
        %defind the plane that is perpenticular to the line defined by the
        %two mid points and passing the post_mid_point
        slope = pre_mid_point-post_mid_point;
        d = dot (slope, pre_mid_point);
        plane_para = cat (2, slope, -d);
        
        %find the point which on the plane and on the cortex and with the
        %longest dis with the pre_mid_point
        allP2Plane_dis = P2Plane_dis (cur_area, plane_para);
        tmp_ind = allP2Plane_dis<4.0;
        P_plane_cortex = cur_area(tmp_ind,:);
        P2post_mp_dis = P2Pset_dis (pre_mid_point, P_plane_cortex);
        [tmp,tmp_ind] = max (P2post_mp_dis);
        P2 = P_plane_cortex(tmp_ind,:);
        %find all the points with the distance to the line defined by the
        %pre_mid_point and the point just found
        P2L_dis = P2Line_dis(cur_area, pre_mid_point, P2);
        cur_division_ind = P2L_dis<step;
        
        %move all the points with the P2line distance lower than the step
        %to the cur_division
        tmp = cur_area (cur_division_ind,:);
        
        if target_area (i) == 2
            tmp(:,1) = -1*tmp(:,1);
        end
        division {i,k} = tmp;
        division_C{i,k} = cur_area_C(cur_division_ind,:);
        plot3(tmp(:,1),tmp(:,2),tmp(:,3),yj_color(k,:),'MarkerSize',10);
        %updating
        pre_mid_point = post_mid_point;
        cur_area (cur_division_ind,:) =[];
        cur_area_C (cur_division_ind,:) =[];
%         size(cur_area,1)
%         pause
    end
    

end

