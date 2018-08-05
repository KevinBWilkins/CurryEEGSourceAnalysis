function [LI_all, cM1, iPM] = plot_act_time_loc (act_ratio_loc_time)
% act_ratio_loc_time (9, 2, 8, 12) = 0.0058;
subjectNo = size (act_ratio_loc_time,1);
for i=1:2
%     ef=squeeze(act_ratio_loc_time(:,i,:,:));
    ef=squeeze(act_ratio_loc_time([1:3,5:9],i,:,:));

    figure
    mean_ef_loc_time=mean(ef,1);
    subplot(1,3,1); 
    tmp = (squeeze(mean_ef_loc_time(:,1:4,:))'); 
    ste = std(ef)./sqrt(subjectNo-1);
    ste_tmp=squeeze(ste(:,1:4,:))';
    errorbar(tmp,ste_tmp)
    cM1=tmp(:,2);
    
    axis ([1 13 0 0.025])
    set(gca,'XTickLabel',-0.6:0.1:-0.1);
    xlabel('Time (s)');
    ylabel ('Activity Ratio');
    legend('S1','M1','PM','SMA')
    grid on
    title ('Contralateral');

    ef_mean = squeeze(mean_ef_loc_time);
    ef_mean = cat (1, ef_mean, ef_mean(end,:));
    subplot(1,3,2)
    surf(squeeze(ef_mean))
    view (0,90)
    axis ([1 13 1 9])
    set(gca,'YTickLabel',{'cS1','cM1','cPM','cSMA','iSMA','iPM','iM1','iS1'})
    set(gca,'XTickLabel',-0.6:0.1:-0.1);
    xlabel('Time (s)');
    if i==1
        title ('SABD')
    elseif i==2
        title ('EF');
    end
    Caxis ([0 0.02])
    colorbar


    subplot(1,3,3)
    tmp1=squeeze(mean_ef_loc_time(:,5:8,:))';
    tmp1 = fliplr(tmp1);
    
    ste_tmp=squeeze(ste(:,5:8,:))';
    ste_tmp = fliplr(ste_tmp);
    errorbar(tmp1,ste_tmp)
    iPM = tmp1 (:,2);
    
    axis ([1 13 0 0.025])
    legend('S1','M1','PM','SMA')
    grid on
    axis ([1 13 0 0.025])
    set(gca,'XTickLabel',-0.6:0.1:-0.1);
    xlabel('Time (s)');
    ylabel ('Activity Ratio');
    title ('Ipsilateral')
    
    
    
    figure (103)
    hold on
   
    contral_all_sub = squeeze(act_ratio_loc_time([1:3,5:9],i,1:4,:));
    ipsi_all_sub = squeeze(act_ratio_loc_time([1:3,5:9],i,5:8,:));
%     contral_all_sub = squeeze(act_ratio_loc_time(:,i,1:4,:));
%     ipsi_all_sub = squeeze(act_ratio_loc_time(:,i,5:8,:));   
    Contral = squeeze (sum (contral_all_sub, 2));
    Ipsi = squeeze (sum (ipsi_all_sub, 2));
    LI = (Contral-Ipsi)./(Contral+Ipsi);
    LI_all (i,:,:) = LI;
    mean_LI= squeeze (mean(LI,1));
    ste_LI = squeeze (std(LI)./sqrt(subjectNo-1));
    subplot (1,2,i)
    errorbar(mean_LI,ste_LI, 'b');  
    
    axis ([1 13 -0.3 0.3])
    set(gca,'XTickLabel',-0.6:0.1:-0.1);
    xlabel('Time (ms)');
    ylabel ('Laterality Index');
    if i ==1 
        title ('SABD');
    else
        title ('EF');
    end
    grid on

end

