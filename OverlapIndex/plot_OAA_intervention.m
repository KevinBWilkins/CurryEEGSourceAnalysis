function  plot_OAA_intervention (pre_OAA, post_OAA)
% act_ratio_loc_time (9, 2, 8, 12) = 0.0058;
% act_ratio_loc_time: subject, task, location, time
subjectNo = size (pre_OAA,1);
figure (101)
hold on
for i_sub = 1: subjectNo
    
        OAA_pre = pre_OAA(i_sub,:);
        OAA_post = post_OAA(i_sub,:);
       
       
        subplot (3,2,i_sub);
        plot([OAA_pre;OAA_post]','.-');  
        legend ('pre-intervention','post-intervention');

        axis ([1 13 0.08 0.25])
        set(gca,'XTickLabel',-0.6:0.1:-0.1);
        xlabel('Time (ms)');
        ylabel ('OAA');
        grid on

    end
end

