% This function calculate the p value of each time between control and
% stroke group. x_control is a 2D array with nxm, n is the number of
% control subjects, and m is the observation number in time domain.
function p = ttest_time (x,y)

[x_no,time_no] = size(x);
[y_no,time_no] = size(y);

for i_time = 1:time_no
    x_cur = x (:,i_time);
    y_cur = y (:,i_time);
    [h, pp] = anova1 ( cat(1,x_cur, y_cur),cat(1,ones(size(x_cur)),zeros(size(y_cur)) ) ,'off');
    p(i_time) = pp{2,6};
    
end