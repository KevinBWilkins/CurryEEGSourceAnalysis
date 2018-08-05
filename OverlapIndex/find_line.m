function cur_line = find_line(seed_ind, cur_area)


x_distance = abs(cur_area-repmat(cur_area(seed_ind,:),size(cur_area,1),1));
x_distance = x_distance (:,1);
line_ind = find (x_distance< 3.8);
cur_line = cur_area (line_ind,:,:);
