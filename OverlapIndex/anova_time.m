% This function calculate the p value of each time between control and
% stroke group. x_control is a 2D array with nxm, n is the number of
% control subjects, and m is the observation number in time domain.
function p = anova_time (x_control, x_stroke)

[control_no,time_no] = size(x_control);
[stroke_no,time_no] = size(x_stroke);

for i_time = 1:time_no
    x = cat(1, x_stroke(:,i_time),x_control(:,i_time));
    grp = cat (1, ones (stroke_no, 1), ones (control_no,1)+1);
    p (i_time) = anovan (x, {grp}, 'display', 'off');
    
end