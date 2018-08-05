function p=statis_time_loc_intervention (pre_loc_time,post_loc_time)

[num_sub, num_task, num_region, num_time ] = size (pre_loc_time);
% [s_num_sub ] = size (act_loc_time_stroke,1);
figure
for i_task = 1: num_task
    for i_region = 1:num_region
        x_pre = squeeze (pre_loc_time(:,i_task,i_region,:));
        x_post = squeeze (post_loc_time(:,i_task,i_region,:));
        p_tmp = ranksum_time (x_pre, x_post);
%         p_tmp = ttest_time (x_pre, x_post);
        p (i_task, i_region,:) = p_tmp;
        p_bin_tmp = zeros (size(p_tmp));
        p_bin_tmp (find (p_tmp<0.05)) = 1;
        p_bin_tmp (find (p_tmp>=0.05 & p_tmp<0.1)) = 0.5;
%         p_bin_tmp (p_tmp>=0.1) = 0;
        p_bin (i_task, i_region,:) = p_bin_tmp;

    end
    subplot(1,2,i_task)
    p_task = squeeze (p_bin(i_task,:,:)); 
    p_task = cat (1, p_task, p_task(end,:));
    surf(p_task)
    view (0,90)
    axis ([1 13 1 9])
    set(gca,'YTickLabel',{'cS1','cM1','cPM','cSMA','iSMA','iPM','iM1','iS1'})
    set(gca,'XTickLabel',-0.6:0.1:-0.1);

%     colorbar
    if i_task ==1
        title ('SABD');
    else
        title ('EF');
    end

end