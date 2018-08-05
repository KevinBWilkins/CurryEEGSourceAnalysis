function  [LI_pre,LI_post]= plot_act_time_loc_intervention (pre_ratio_loc_time, post_ratio_loc_time)
% act_ratio_loc_time (9, 2, 8, 12) = 0.0058;
% act_ratio_loc_time: subject, task, location, time
subjectNo = size (pre_ratio_loc_time,1);
for i_sub = 1: subjectNo
    for i_task=1:2
        
        pre=squeeze(pre_ratio_loc_time(i_sub,i_task,:,:));
        post=squeeze(post_ratio_loc_time(i_sub,i_task,:,:));
        pre = cat (1, pre, pre(end,:));
        post = cat (1, post, post(end,:));
        figure (100)
        hold on 

        subplot(3,8,(i_sub-1)*4+(i_task-1)*2+1); 
%         disp (num2str((i_sub-1)*4+(i_task-1)*2+1))
        surf(pre)
        view (0,90)
        axis ([1 13 1 9])
        set(gca,'YTickLabel',{'cS1','cM1','cPM','cSMA','iSMA','iPM','iM1','iS1'})
        set(gca,'XTickLabel',-0.6:0.1:-0.1);
        colorbar

        subplot(3,8,(i_sub-1)*4+(i_task-1)*2+2); 
        surf(post)
%         disp (num2str((i_sub-1)*4+(i_task-1)*2+2))
        view (0,90)
        axis ([1 13 1 9])
        set(gca,'YTickLabel',{'cS1','cM1','cPM','cSMA','iSMA','iPM','iM1','iS1'})
        set(gca,'XTickLabel',-0.6:0.1:-0.1);
        colorbar

        figure (101)
        hold on
        contral_pre_sub = sum(pre(1:4,:), 1);
        ipsi_pre_sub = sum(pre(5:8,:), 1);   
        contral_post_sub = sum(post(1:4,:), 1);
        ipsi_post_sub = sum(post(5:8,:), 1);   
        
        
        LI_pre ( (i_task-1)* subjectNo + i_sub, : ) = (contral_pre_sub-ipsi_pre_sub)./(contral_pre_sub+ipsi_pre_sub);
        LI_post ( (i_task-1)* subjectNo + i_sub, : ) = (contral_post_sub-ipsi_post_sub)./(contral_post_sub+ipsi_post_sub);
       
        subplot (3,4,(i_sub-1)*2+i_task)
        plot([LI_pre;LI_post]','.-');  
        legend ('pre-intervention','post-intervention');

        axis ([1 13 -0.7 0.7])
        set(gca,'XTickLabel',-0.6:0.1:-0.1);
        xlabel('Time (ms)');
        ylabel ('Laterality Index');
        if i_task ==1 
            title ('SABD');
        else
            title ('EF');
        end
        grid on

    end
end

