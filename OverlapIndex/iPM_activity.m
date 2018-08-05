function iPM=iPM_activity (act_time_loc)

[sub_num, task_num, loc_num, time_num]=size(act_time_loc);
for i = 1: sub_num
    for j=1:task_num
        iPM (i,j) = max( squeeze(act_time_loc(i,j,7,:)) );
    end
end