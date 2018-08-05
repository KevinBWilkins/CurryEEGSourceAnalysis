% This function calculate the p value of each time between control and
% stroke group. x_control is a 2D array with nxm, n is the number of
% control subjects, and m is the observation number in time domain.
function p = ranksum_time (x,y)

[x_no,time_no] = size(x);
[y_no,time_no] = size(y);

for i_time = 1:time_no
    x_cur = x (:,i_time);
    y_cur = y (:,i_time);
    [ pp,h] = ranksum (x_cur, y_cur);
    p(i_time) = pp;
    
end