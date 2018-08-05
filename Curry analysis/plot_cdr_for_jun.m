

%20. Plot cortex activity within time bins of 50 ms
if areanum==7

    %20a. Plot
    figure(20+4*areanum)
    set(20+4*areanum,'Color',[1 1 1])
    figure(21+4*areanum)
    set(21+4*areanum,'Color',[1 1 1])

    index=0;
    %8 subplots

    %if starts at -2000, then add 8*timeinc
    %else starts at -1000, add 4*timeinc
    if strfind(cur_task,'mri')
        addtime1 = 4*timeinc+1;
        addtime2 = -4*timeinc;
    else
        addtime1 = 8*timeinc;
        addtime2 = 0*timeinc;;
    end

    for timenow = timestart+addtime1:timeinc:timestop+addtime2
        index=index+1;

        %8 or 16 subplots
        figure(20+4*areanum)
        %subplot(4,4,index)
        subplot(2,4,index)
        imshow(imA_resize2,'XData',[mincortX maxcortX],'YData',[maxcortY mincortY])
        axis([mincortX maxcortX mincortY maxcortY])
        axis off
        hold on

        timestart_seg = timenow;
        timestop_seg = timenow + timeinc-1;
        if timestop_seg > timestop
            timestop_seg = timestop;
        end

        norm_mean_loc = mean(norm_mean_locs(:,timestart_seg:timestop_seg),2);


        %SURFACE PLOT
        C = mean(newtargetC(:,timestart_seg:timestop_seg),2);
        CI = griddata(newtargetL(:,1),newtargetL(:,2),C,newrot_XI,newrot_YI);
        figure(5)
        h3 = surface('XData',newrot_XI,'YData',newrot_YI,'ZData',newrot_ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6);
        vaxis2 = caxis;
        close(5)

        figure(20+4*areanum)
        subplot(2,4,index)
        h3 = surface('XData',newrot_XI,'YData',newrot_YI,'ZData',newrot_ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6);

        set(h3,'FaceAlpha','flat','AlphaData',new_mask_perim_filled_resize);
        set(gca,'ALim',[0 1]);


        view(0,90)
        TR = timeRange_total{i};
        title([num2str(TR(1)+(TR(2)-TR(1))*(timestart_seg-1)/(timePoints_total(i)-1),5),' to ',num2str(TR(1)+(TR(2)-TR(1))*(timestop_seg-1)/(timePoints_total(i)-1),5),'ms'],'Fontsize',14,'FontWeight','bold')
        axis off
        caxis([0 1])
        %caxis(vaxis2)
        axis equal
        axis xy
        axis tight

        a1 = gca;
        a2 = axes;
        caxis([0 1])
        %caxis(vaxis2);
        colorbar
        set(a2,'Position',get(a1,'Position'));
        axis off;

        
    end
    figure(20+4*areanum)
    hgsave(gcf,[pwd,'\figures\',subjectName,'\',subjectName,'_',cur_task,'_8brains.fig']);
end