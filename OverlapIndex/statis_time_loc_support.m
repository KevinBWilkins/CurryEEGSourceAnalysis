function p=statis_time_loc_support (on_loc_time,off_loc_time)

[num_time,num_freq,num_sub, num_task] = size (on_loc_time);
% [s_num_sub ] = size (act_loc_time_stroke,1);
figure
for i_task = 1:num_task
    for i_freq = 1:num_freq
        x_on = squeeze (on_loc_time(:,i_freq,:,i_task));
        x_off = squeeze (off_loc_time(:,i_freq,:,i_task));
        p_tmp = ttest_time (x_on, x_off);
        p (i_task, i_freq,:) = p_tmp;
        p_bin_tmp = zeros (size(p_tmp));
        p_bin_tmp (find (p_tmp<0.05)) = 1;
        p_bin_tmp (find (p_tmp>=0.05 & p_tmp<0.1)) = 0.5;
    %         p_bin_tmp (p_tmp>=0.1) = 0;
        p_bin (i_task, i_freq,:) = p_bin_tmp;

    end
    subplot(1,2,i_task)
    p_task = squeeze (p_bin(i_task,:,:)); 
    p_task = cat (1, p_task, p_task(end,:));
    surf(p_task)
    view (0,90)
    axis ([1 13 1 72])

    %     colorbar
    if i_task ==1
        title ('close');
    else
        title ('open');
    end
end

