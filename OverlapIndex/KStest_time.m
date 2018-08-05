% This function calculate the p value of each time between control and
% stroke group. x_control is a 2D array with nxm, n is the number of
% control subjects, and m is the observation number in time domain.
function h = KStest_time (x,y)

[x_no,time_no] = size(x);
[y_no,time_no] = size(y);

for i_time = 1:time_no
    x_cur = x (:,i_time);
    y_cur = y (:,i_time);
    [hx] = kstest (x_cur,[],0.1);
    [hy] = kstest (y_cur,[],0.1);
    h(i_time) = hx*hy;
    
end