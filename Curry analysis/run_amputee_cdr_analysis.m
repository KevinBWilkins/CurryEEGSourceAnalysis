function [results] = run_amputee_cdr_analysis()



%load all amputee files

reject_electrodes = [];
%--------------------------------------------------------------------------
% %LS
% %amputeesubjectnames = ['LS1','LS2','LSpost1','LSpost2'];
% 
% %tasknames = ['abd','hcl'];
% 
% amputeepregrandlist = {'LSpost1_abd_un_mri3','LSpost1_ef_un_mri2','LSpost1_hcl_un_mri','LSpost1_ho_un_mri',...
%     };
% amputeepregrandlist2 = {'LS1_abd2','LS1_ef2','LS1_hcl2','LS2_ho2',...
%     };
% 
% amputeepostgrandlist = {'LSpost1_abd_un_mri3','LSpost1_ef_un_mri2','LSpost1_hcl_un_mri','LSpost1_ho_un_mri',...
%     };
% amputeepostgrandlist2 = {'LSpost2_abd_mri','LSpost2_ef_mri','LSpost2_hcl_mri','LSpost2_ho_mri',...
%     };
% 
% sides = {'left','right'};       %intact = first, amputated = second
% %sides = {'both','both'};       %intact = first, amputated = second
% reject_electrodes = [14 0;14 0;15 16;16 0];
%--------------------------------------------------------------------------
%AK
%amputeesubjectnames = ['AKpre1','AKpre2','AKpost1','AKpost2','LS1','LS2','LSpost1','LSpost2'];

%tasknames = ['abd_new','ee_new','hcl_new'];

% amputeepregrandlist = {'AKpre1_abd_un_new','AKpre1_ee_un_new2','AKpre1_hcl_un_new2',...
%     };
amputeepregrandlist = {'AKpost2_abd_un_new2','AKpost2_ee_un_new','AKpost2_hcl_un_new',...
    };
amputeepregrandlist2 = {'AKpre2_abd_new2','AKpre2_ee_new','AKpre2_hcl_new',...
    };

amputeepostgrandlist = {'AKpost2_abd_un_new2','AKpost2_ee_un_new','AKpost2_hcl_un_new'...
    };
amputeepostgrandlist2 = {'AKpost1_abd_new3','AKpost1_ee_new','AKpost1_hcl_new2'...
    };
sides = {'right','left'};       %intact = first, amputated = second
reject_electrodes = [0 0 0 0 0;0 0 0 0 16;0 0 0 0 0;12 13 14 15 16];
%--------------------------------------------------------------------------
% %CM
% %amputeesubjectnames = ['CM103006','CM110206'];
%
% %tasknames = ['abd','hcl'];
%
% amputeepregrandlist = {...
%     };
% amputeepregrandlist2 = {...
%     };
%
% amputeepostgrandlist = {'CM103006_abd_un','CM103006_hcl_un'...
%     };
% amputeepostgrandlist2 = {'CM110206_abd','CM110206_hcl'...
%     };
%--------------------------------------------------------------------------
close all;

dmean_m = [];
dmean_m2 = [];
init_max_dist_btwn_clusters = 15;

%iterate through different numbers of clusters (k) 
k_m = [2:2:2];

for i_k = 1:length(k_m)
    for i_iter = 1:1

        %loop through all lists
        for i_list = 1:4

            %initialization
            grandlist = [];

            if i_list==1        %PRE-TR INTACT SIDE
                grandlist = amputeepregrandlist;
                plot_title = 'pre-TR int';
                fig_num = 1;
            elseif i_list==2    %POST-TR INTACT SIDE
                grandlist = amputeepostgrandlist;
                plot_title = 'post-TR int';
                fig_num = 1;
            elseif i_list==3    %PRE-TR AMPUTATED SIDE
                grandlist = amputeepregrandlist2;
                plot_title = 'pre-TR amp';
                fig_num = 2;
            elseif i_list==4    %POST-TR AMPUTATED SIDE
                grandlist = amputeepostgrandlist2;
                plot_title = 'post-TR amp';
                fig_num = 2;
            end



            %go through tasks within each list
            for i_task=1:length(grandlist)


                %read in file
                load(['C:\Albert Chen\matlab\Curry analysis\cdr\',grandlist{i_task},'_cdr.mat']);

                if strfind(grandlist{i_task},'abd')
                    m_color = 'b';
                    task_name = 'abd';
                elseif strfind(grandlist{i_task},'ee')
                    %m_color = 'g';
                    m_color = 'b';
                    task_name = 'ee';
                elseif strfind(grandlist{i_task},'ef')
                    %m_color = 'g';
                    m_color = 'b';
                    task_name = 'ef';
                elseif strfind(grandlist{i_task},'hcl')
                    %m_color = 'm';
                    m_color = 'b';
                    task_name = 'hcl';
                elseif strfind(grandlist{i_task},'ho')
                    %m_color = 'c';
                    m_color = 'b';
                    task_name = 'ho';
                end

                %recenter image by calculating x image shift
                % x_image_shift = mean(cdr.newtargetL(:,1));
                % x_image_shift = (cdr.maxcortX-cdr.mincortX)/2 * mean(cdr.norm_loc_cortex2(:,1)) + (cdr.mincortX+cdr.maxcortX)/2;
                x_image_shift = (cdr.mincortX+cdr.maxcortX)/2;

                cdr.newtargetL = cdr.newtargetL-x_image_shift;
                cdr.rot_new_perim = cdr.rot_new_perim - x_image_shift;
                
                cdr.mincortX = cdr.mincortX-x_image_shift;
                cdr.maxcortX = cdr.maxcortX-x_image_shift;

                %-----------------------
                %display cortex image
                figure(fig_num+20*(i_k-1))
                set(gcf,'Color',[1 1 1]);

                subplot(2,length(grandlist),mod(i_list-1,2)*length(grandlist)+i_task)
                imshow(cdr.imA_resize2,'XData',[cdr.mincortX cdr.maxcortX],'YData',[cdr.maxcortY cdr.mincortY])
                hold on;

                %-----------------------
                %display activesources
                threshold = cdr.newThreshold;
                activesources = double(cdr.activesources>0);

                %use same threshold for all subjects based on max activation
%                             activesources = double(cdr.targetCI>(.7));
%                             threshold = .7;


                
                %----------------------------------------------------
                %plot single COG over single hemisphere

                if isfield(cdr,'left_group')
                    left_targetCI = cdr.targetCI(1:sum(cdr.left_group>0),:);
                    right_targetCI = cdr.targetCI((sum(cdr.left_group>0)+1):end,:);
                    left_newtargetL = cdr.newtargetL(1:sum(cdr.left_group>0),:);
                    right_newtargetL = cdr.newtargetL((sum(cdr.left_group>0)+1):end,:);

                    left_activesources = activesources(1:sum(cdr.left_group>0),:);
                    right_activesources = activesources((sum(cdr.left_group>0)+1):end,:);

                    if isequal(lower(sides{fig_num}),'right')
                        targetCI = left_targetCI;
                        activesources2 = left_activesources;
                        newtargetL = left_newtargetL;
                    else
                        targetCI = right_targetCI;
                        activesources2 = right_activesources;
                        newtargetL = right_newtargetL;
                    end
                else
                    if isequal(sides{fig_num},'right')
                        %find sources on contralateral, left side of brain
                        group = (cdr.newtargetL(:,1)<0).*(1:length(cdr.newtargetL(:,1)))';
                    elseif isequal(sides{fig_num},'left')
                        %find sources on contralateral, right side of brain
                        group = (cdr.newtargetL(:,1)>0).*(1:length(cdr.newtargetL(:,1)))';
                    else
                        %use both sides
                        group = 1:length(cdr.targetCI(:,1));
                    end
                    group = group(group>0);
                    targetCI = cdr.targetCI(group,:);
                    activesources2 = activesources(group,:);
                    newtargetL = cdr.newtargetL(group,:);
                end

                totalstrengths = sum(targetCI.*activesources2,1);
                for coord = 1:3
                    [t_m,t_n] = size(targetCI);
                    targetL_M = repmat(newtargetL(:,coord),1,t_n);
                    weighted_locs = targetL_M.*activesources2.*targetCI;
                    mean_locs(coord,1:t_n) = sum(weighted_locs,1)./totalstrengths;
                end

                mean_locs_window = mean_locs(:,(cdr.plotstop-18):(cdr.plotstop+8));
                %remove NaN's if any
                num_nan = sum(isnan(mean_locs_window),2);
                mean_locs_window(isnan(mean_locs_window)) = 0;
                single_mean_loc = mean(mean_locs_window,2)';
                %correct for mean
                single_mean_loc = single_mean_loc*length(mean_locs_window(1,:))/(length(mean_locs_window(1,:)) - num_nan(1));

                %plot single COG
                if i_list == 3
                    plot(single_mean_loc(1),single_mean_loc(2),'Markersize',30,'color','r','LineStyle','.','Linewidth',3)
                elseif i_list == 4
                    plot(single_mean_loc(1),single_mean_loc(2),'Markersize',30,'color','b','LineStyle','.','Linewidth',3)
                else
                    plot(single_mean_loc(1),single_mean_loc(2),'Markersize',30,'color','c','LineStyle','.','Linewidth',3)
                end

                single_mean_loc_m{i_list,i_task} = single_mean_loc;
                mirror_single_mean_loc = single_mean_loc;
                mirror_single_mean_loc(1) = -mirror_single_mean_loc(1);
                mirror_single_mean_loc_m{i_list,i_task} = mirror_single_mean_loc;

                %plot mirror COG for intact side if plotting amputated side
                if i_list>2
                    mirror_intact_mean_loc = mirror_single_mean_loc_m{i_list-2,i_task};
                    plot(mirror_intact_mean_loc(1),mirror_intact_mean_loc(2),'Markersize',30,'color','c','LineStyle','.','Linewidth',3)
                else
                    plot(mirror_single_mean_loc(1),mirror_single_mean_loc(2),'Markersize',30,'color','c','LineStyle','.','Linewidth',3)
                end

                %-----------------------
                %quantify number of active sources (size of activity) before and after TR
                numsources = sum(activesources2,1);
                totalnumsources = length(newtargetL(:,1));
                
                activesources_70 = double(targetCI>(.7));
                numsources_70 = sum(activesources_70);
                                
                figure(5)
                set(gcf,'Color',[1 1 1]);
                %title('Sizes of activity')
                %subplot(2,length(grandlist),mod(i_list-1,2)*length(grandlist)+i_task)
                subplot(1,length(grandlist),i_task)
                hold on;
                if fig_num == 1
                    plot(-500:(800/(length(numsources)-1)):300,numsources/totalnumsources,'c','Linewidth',2);
                    %plot(-500:(800/(length(numsources_70)-1)):300,numsources_70,'b','Linewidth',2,'Linestyle','--');
                else
                    if i_list == 3
                        plot(-500:(800/(length(numsources)-1)):300,numsources/totalnumsources,'r','Linewidth',2);
                    elseif i_list == 4
                        plot(-500:(800/(length(numsources)-1)):300,numsources/totalnumsources,'b','Linewidth',2);
                    end
                    %plot(-500:(800/(length(numsources_70)-1)):300,numsources_70,'r','Linewidth',2,'Linestyle','--');
                end
                set(gca,'XLim',[-500 300]);
                xlabel('time (ms)');
                ylabel('size/(total size of ROI)');
                mean_numsources(i_list,i_task) = mean(numsources)/totalnumsources;
                mean_numsources_70(i_list,i_task) = mean(numsources_70)/totalnumsources;
                
                %-----------------------
                %quantify mean activity of active sources (strength of activity) before and after TR
                meanstrengths = mean(targetCI,1);
                
                meanstrengths_70 = mean(targetCI.*double(targetCI>(.7)),1);
                                
                figure(6)
                set(gcf,'Color',[1 1 1]);
                %title('Mean strengths of activity')
                %subplot(2,length(grandlist),mod(i_list-1,2)*length(grandlist)+i_task)
                subplot(1,length(grandlist),i_task)
                hold on;
                if fig_num == 1
                    plot(-500:(800/(length(meanstrengths)-1)):300,meanstrengths,'c','Linewidth',2);
                    %plot(-500:(800/(length(meanstrengths_70)-1)):300,meanstrengths_70,'b','Linewidth',2,'Linestyle','--');
                else
                    if i_list == 3
                        plot(-500:(800/(length(meanstrengths)-1)):300,meanstrengths,'r','Linewidth',2);
                    elseif i_list == 4
                        plot(-500:(800/(length(meanstrengths)-1)):300,meanstrengths,'b','Linewidth',2);
                    end
                    %plot(-500:(800/(length(meanstrengths_70)-1)):300,meanstrengths_70,'r','Linewidth',2,'Linestyle','--');
                end
                set(gca,'XLim',[-500 300]);

                mean_meanstrengths(i_list,i_task) = mean(meanstrengths);
                mean_meanstrengths_70(i_list,i_task) = mean(meanstrengths_70);
                
                xlabel('time (ms)');
                ylabel('strength');
                %----------------------------------------------------------

                %repeat analysis if the shortest distance between clusters for different days is > 10mm
                %stop analysis if repeat more than 10 times
                distancetoobig = 1;
                timesrepeated = 0;
                timesplotfig = 0;
                max_dist_btwn_clusters = init_max_dist_btwn_clusters;
                while distancetoobig && timesrepeated < 10
                    timesplotfig = timesplotfig + 1;
                    %----------------------------------------------------
                    %use weighted k-means algorithm (wkmeans.m) to get multiple clusters and COG's
                    %X = locations, k = number of clusters, W = weights

                    %   IDX = WKMEANS(X, K, W) partitions the points in the N-by-P data matrix
                    %   X with an N-by-1 weight matrix W into K clusters.  This partition minimizes the sum, over all
                    %   clusters, of the weighted within-cluster sums of point-to-cluster-centroid
                    %   distances.  Rows of X correspond to points, columns correspond to
                    %   variables.  WKMEANS returns an N-by-1 vector IDX containing the
                    %   cluster indices of each point.  By default, WKMEANS uses weighted squared
                    %   Euclidean distances.
                    %
                    %   Modified from KMEANS so that the weights in W are used in the distance
                    %   calculation when finding the centroids.
                    %
                    %   Example:
                    %
                    %       X = [randn(20,2)+ones(20,2); randn(20,2)-ones(20,2)];
                    %       W = [randn(20,1); randn(20,1)];
                    %       [cidx, ctrs] = wkmeans(X, 2, W);
                    %       plot(X(cidx==1,1),X(cidx==1,2),'r.', ...
                    %            X(cidx==2,1),X(cidx==2,2),'b.', ctrs(:,1),ctrs(:,2),'kx');

                    %make matrices of active sources and corresponding weights (strengths)
                    indices = repmat((1:length(newtargetL(:,1)))',1,27).*activesources2(:,(cdr.plotstop-8):(cdr.plotstop+18));
                    indices = sum(indices,2);
                    indices2 = (1:length(newtargetL(:,1)))';
                    indices2 = indices2(indices>0);
                    X = newtargetL(indices2,1:2);
                    W = targetCI(indices2,(cdr.plotstop-8):(cdr.plotstop+18));
                    W = mean(W,2);

                    k = min(k_m(i_k),length(X(:,1)));
                    [cidx, ctrs] = wkmeans(X, k, W);

                    
                    %save active sources
                    indices_m{i_list,i_task} = indices2;
                    newtargetL_m{i_list,i_task} = newtargetL;
                    newtargetL_m2{i_list,i_task} = cdr.newtargetL;
                    W_m{i_list,i_task} = W;
                    
                    %----------------------------------------------------
                    indices3 = repmat((1:length(cdr.newtargetL(:,1)))',1,27).*activesources(:,(cdr.plotstop-8):(cdr.plotstop+18));
                    indices3 = sum(indices3,2);
                    indices4 = (1:length(cdr.newtargetL(:,1)))';
                    indices4 = indices4(indices3>0);

                    %Plot active locations- either in one color or as surface with strengths
                    figure(fig_num+20*(i_k-1))
                    subplot(2,length(grandlist),mod(i_list-1,2)*length(grandlist)+i_task)
%                     %plot active sources and COG's of the clusters
                    %plot3(newtargetL(indices2,1),newtargetL(indices2,2),newtargetL(indices2,3),'Color',m_color,'Markersize',10,'LineStyle','.');
                    %plot3(Subject.mean_locs(1,:),Subject.mean_locs(2,:),cd
                    %r.newmaxrot_Z*ones(1,length(Subject.mean_locs(1,:))),'Color',m_color,'Markersize',7,'LineStyle','.')
                    
                    
                    %
                    plot3(cdr.newtargetL(indices4,1),cdr.newtargetL(indices4,2),cdr.newtargetL(indices4,3),'Color','m','Markersize',10,'LineStyle','.');
                    

%                     if timesplotfig == 1
%                         newrot_XI = cdr.norm_newrot_XI.*(cdr.maxcortX-cdr.mincortX)/2 + (cdr.mincortX+cdr.maxcortX)/2;
%                         newrot_YI = cdr.norm_newrot_YI.*(cdr.maxcortY-cdr.mincortY)/2 + (cdr.mincortY+cdr.maxcortY)/2;
%                         newrot_ZI = cdr.newrot_ZI;
%                         C = mean(cdr.targetCI(:,(cdr.plotstop-18):(cdr.plotstop+8)),2);
%                         CI = griddata(cdr.newtargetL(:,1),cdr.newtargetL(:,2),C,newrot_XI,newrot_YI);
% 
%                         %plot dummy figure
%                         figure(100)
%                         h3 = surface('XData',newrot_XI,'YData',newrot_YI,'ZData',newrot_ZI,'CData',CI,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6);
%                         v_caxis = caxis;
%                         close(100);
% 
%                         figure(fig_num+20*(i_k-1))
%                         subplot(2,length(grandlist),mod(i_list-1,2)*length(grandlist)+i_task)
%                         h3 = surface('XData',newrot_XI,'YData',newrot_YI,'ZData',newrot_ZI,'CData',CI,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.6);
%                         set(h3,'FaceAlpha','flat','AlphaData',cdr.new_mask_perim_filled_resize);
%                         set(gca,'ALim',[0 1]);
%                         caxis(v_caxis);
%                     end
                    
                    %plot(ctrs(:,1),ctrs(:,2),'kx','Markersize',10,'Linewidth',3);
                    % X(cidx==1,1),X(cidx==1,2),'r.', X(cidx==2,1),X(cidx==2,2),'b.',


                    %plot single mirrored location for cluster with most members
                    mirror_ctrs = ctrs;
                    mirror_ctrs(:,1) = -mirror_ctrs(:,1);

                    if strfind(grandlist{i_task},'abd')
                        %use cluster closest to midline
                        [y_min i_min] = min(abs(mirror_ctrs(:,1)));
                        single_mirror_ctr = mirror_ctrs(i_min,:);
                    elseif strfind(grandlist{i_task},'hcl') & (i_list < 3) & strfind(grandlist{i_task},'AK') & k < 4
                        %use cluster furthest from midline
                        [y_max i_max] = max(abs(mirror_ctrs(:,1)));
                        single_mirror_ctr = mirror_ctrs(i_max,:);
                    elseif strfind(grandlist{i_task},'ho') & (i_list < 3) & strfind(grandlist{i_task},'LS') & k < 3
                        %use cluster furthest from midline
                        [y_max i_max] = max(abs(mirror_ctrs(:,1)));
                        single_mirror_ctr = mirror_ctrs(i_max,:);
                    else

                        %use mirror_ctr cluster with most members
                        unique_members = unique(cidx);
                        for i_uni = 1:length(unique_members)
                            num_members(i_uni) = length(find(cidx==unique_members(i_uni)));
                        end
                        [y_max i_max] = max(num_members);
                        single_mirror_ctr = mirror_ctrs(i_max,:);
                    end

                    mirror_ctrs_m{i_list,i_task} = single_mirror_ctr;
                    if i_list < 3
                        %plot(single_mirror_ctr(1),single_mirror_ctr(2),'Markersize',30,'color','c','LineStyle','.','Linewidth',3)
                    else
                        mirror_intact_single_cluster = mirror_ctrs_m{i_list-2,i_task};
                        %plot(mirror_intact_single_cluster(1),mirror_intact_single_cluster(2),'Markersize',30,'color','c','LineStyle','.','Linewidth',3)
                    end
                    %----------------------------------------------------------

                    %fix view
                    view(0,90)
                    axis on

                    %create grid for show purposes
                    %     set(gca,'XTick',30:5:80)
                    %     set(gca,'YTick',70:5:160)
                    %grid on
                    axis xy


                    axis equal
                    axis([cdr.mincortX cdr.maxcortX cdr.mincortY cdr.maxcortY])
                    axis tight
                    %     axis([30 80 70 160])
                    title([plot_title,'- ',task_name]);


                    %save information
                    cidx_matrix{i_list,i_task} = cidx;
                    ctrs_matrix{i_list,i_task} = ctrs;
                    threshold_m{i_list,i_task} = threshold;
                    


                    %------------------------------------------------------------------
                    %Figure out distances between clusters
                    %limit differences to 10mm for multiple COG
                    %do not limit for single COG

                    %multiple COG analysis
                    %         i2 = 0;
                    %         dmean = [];
                    %         for i = 1:2:3
                    %             i2 = i2+1;
                    %             for j=1:i_task
                    %
                    %                 %find distances between two lists of centers
                    %                 [d{i2,j}, xshift{i2,j}, yshift{i2,j}, final_ind_m{i2,j}] = distlistfun(ctrs_matrix{i,j},ctrs_matrix{i+1,j},1);
                    %
                    %                 %mean distance for each task
                    %                 dmean(i2,j) = mean(d{i2,j});
                    %                 xshiftmean(i2,j) = mean(xshift{i2,j});
                    %                 yshiftmean(i2,j) = mean(yshift{i2,j});
                    %             end
                    %         end

                    %COG analysis
                    if i_list > 2
                        i2 = i_list-2;


                        %MULTIPLE COG- BETWEEN AMPUTATED AND SINGLE MIRRORED INTACT CLUSTER
                        %find distances between two lists of centers

                        [d{i2,i_task}, xshift{i2,i_task}, yshift{i2,i_task}, final_ind_m{i2,i_task}] = distlistfun(mirror_ctrs_m{i_list-2,i_task},ctrs_matrix{i_list,i_task},max_dist_btwn_clusters);

                        %mean distance for each task
                        dmean(i2,i_task) = mean(d{i2,i_task});
                        xshiftmean(i2,i_task) = mean(xshift{i2,i_task});
                        yshiftmean(i2,i_task) = mean(yshift{i2,i_task});


                        %SINGLE COG- BETWEEN MIRRORED INTACT SIDE AND AMPUTATED
                        %find distances between post mirrored intact side and both amputated sides
                        mirror_intact_mean_loc = single_mean_loc_m{i2,i_task};
                        mirror_intact_mean_loc(1) = -mirror_intact_mean_loc(1);
                        [d_single{i2,i_task}, xshift_single{i2,i_task}, yshift_single{i2,i_task}, final_ind_single_m{i2,i_task}] = distlistfun(mirror_intact_mean_loc,single_mean_loc_m{i_list,i_task},0);

                        %mean distance for each task
                        dmean_single(i2,i_task) = mean(d_single{i2,i_task});
                        xshiftmean_single(i2,i_task) = mean(xshift_single{i2,i_task});
                        yshiftmean_single(i2,i_task) = mean(yshift_single{i2,i_task});
                        
                        
                        %SINGLE COG- BETWEEN AMPUTATED SIDES
                        %find distances between post amputated side and pre amputated side
                        if i_list ==4
                            [d_single_amp{1,i_task}, xshift_single_amp{1,i_task}, yshift_single_amp{1,i_task}, final_ind_single_amp_m{1,i_task}] = distlistfun(single_mean_loc_m{i_list-1,i_task},single_mean_loc_m{i_list,i_task},0);
                        
                            dmean_single_amp(1,i_task) = (d_single_amp{1,i_task});
                        end
                        
                        
                        

                        %check to make sure there are real distances
                        if (sum(sum(isnan(dmean)))) == 0
                            distancetoobig = 0;
                            max_dist_btwn_clusters_m(i_iter) = max_dist_btwn_clusters;
                            %reset max_dist_btwn_clusters to init_max_dist_btwn_clusters;
                            max_dist_btwn_clusters = init_max_dist_btwn_clusters;
                        else
                            fprintf(1,'distance too big between clusters\n');
                            max_dist_btwn_clusters = max_dist_btwn_clusters+5;
                            timesrepeated = timesrepeated+1;
                        end
                    else
                        distancetoobig = 0;
                    end




                    %----------------------------------------------------------
                    %display EMGs

                    if timesplotfig == 1
                        %load EMG .mat file
                        [name, taskname] = strtok(grandlist{i_task},'_');
                        taskname = [name, '_', strtok(taskname,'_')];
                        if i_list <= 2
                            taskname = [taskname,'_un'];
                        end
                        load(['C:\Albert Chen\Torque\',taskname,'_emg.mat']);

                        figure(10+fig_num)
                        set(gcf,'Color',[1 1 1]);

                        subplot(2,length(grandlist),mod(i_list-1,2)*length(grandlist)+i_task)

                        if ~isempty(reject_electrodes)
                            for i_rej = 1:length(reject_electrodes(i_list,:))
                                if reject_electrodes(i_list,i_rej)~=0
                                    emg.whole_maxEMGs(reject_electrodes(i_list,i_rej)) = 0;
                                    emg.maxEMGs(reject_electrodes(i_list,i_rej)) = 0;
                                end
                            end
                        end

                        disp_electrodes = [1:3 5:16];
                        bar(1:15,emg.whole_maxEMGs(disp_electrodes),'k');
                        hold on;

                        bar(1:15,emg.maxEMGs(disp_electrodes),m_color,'barwidth',.5);

                        x = 1:15;
                        y = emg.maxEMGs(disp_electrodes);
                        emgnames = emg.names(disp_electrodes);

                        ylim = get(gca,'YLim');

                        %label individual bars with EMG name
                        for i_text = 1:length(x)
                            if isequal(emgnames{i_text},'FCD')
                                emgnames{i_text} = 'FCU';
                            elseif isequal(emgnames{i_text},'ECD')
                                emgnames{i_text} = 'EDC';
                            elseif isequal(emgnames{i_text},'PMJ_m_hcl')
                                emgnames{i_text} = 'PMJ_mhcl';
                            elseif isequal(emgnames{i_text},'PMJ_m_ho')
                                emgnames{i_text} = 'PMJ_mho';
                            elseif isequal(emgnames{i_text},'PMJ_m_ef')
                                emgnames{i_text} = 'PMJ_mef';
                            elseif isequal(emgnames{i_text},'TRAPupper')
                                 emgnames{i_text} = 'TRAPS';
                            end
                            text(x(i_text),-ylim(2)/25,emgnames{i_text},'Rotation',-90);
                        end

                        set(gca,'XLim',[0 16],'YLim',[0 ylim(2)],'xticklabel',{' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',})
                        title([plot_title,'- ',task_name]);
                        
                        %save EMG info
                        totalEMG_m(i_list,i_task) = sum(emg.maxEMGs(disp_electrodes));
                        
                    end


                end
            end

            if i_list == 4

                if sum(sum(isnan(dmean))) == 0
                    %distance between multiple COG for same side
                    dmean_m{i_iter} = dmean;
                    dmean_m2(i_iter,:) = [dmean(1,1:end) dmean(2,1:end)]

                    xshiftmean_m{i_iter} = xshiftmean;
                    xshiftmean_m2(i_iter,:) = [xshiftmean(1,1:end) xshiftmean(2,1:end)];

                    yshiftmean_m{i_iter} = yshiftmean;
                    yshiftmean_m2(i_iter,:) = [yshiftmean(1,1:end) yshiftmean(2,1:end)];

                    %distance between single COG for amputated and mirrored intact side
                    dmean_single_m{i_iter} = dmean_single;
                    dmean_single_m2(i_iter,:) = [dmean_single(1,1:end) dmean_single(2,1:end)]

                    xshiftmean_single_m{i_iter} = xshiftmean_single;
                    xshiftmean_single_m2(i_iter,:) = [xshiftmean_single(1,1:end) xshiftmean_single(2,1:end)];

                    yshiftmean_single_m{i_iter} = yshiftmean_single;
                    yshiftmean_single_m2(i_iter,:) = [yshiftmean_single(1,1:end) yshiftmean_single(2,1:end)];
                else
                    dmean_m{i_iter} = 0;
                    dmean_m2(i_iter,:) = 0;
                    xshiftmean_m{i_iter} = 0;
                    xshiftmean_m2(i_iter,:) = 0;
                    yshiftmean_m{i_iter} = 0;
                    yshiftmean_m2(i_iter,:) = 0;
                    dmean_single_m{i_iter} = 0;
                    dmean_single_m2(i_iter,:) = 0;
                    xshiftmean_single_m{i_iter} = 0;
                    xshiftmean_single_m2(i_iter,:) = 0;
                    yshiftmean_single_m{i_iter} = 0;
                    yshiftmean_single_m2(i_iter,:) = 0;
                end
            end
        end
    end


    %plot mean sizes of activity
    figure(7)
    mean_numsources
    for i_bar = 1:length(mean_numsources(1,:))
        subplot(1,length(mean_numsources(1,:)),i_bar)
        bar(mean_numsources(:,i_bar),'group');
    end

    %plot mean strengths of activity
    figure(8)
    mean_meanstrengths
    for i_bar = 1:length(mean_meanstrengths(1,:))
        subplot(1,length(mean_meanstrengths(1,:)),i_bar)
        bar(mean_meanstrengths(:,i_bar),'group');
    end
    
    %plot overlap of amputated active sources
    figure(9)
    set(gcf,'Color',[1 1 1]);
    for i_task = 1:length(grandlist)
        
        %subplot(1,length(grandlist),i_task)
        subplot(1,4,i_task)
        
        imshow(cdr.imA_resize2,'XData',[cdr.mincortX cdr.maxcortX],'YData',[cdr.maxcortY cdr.mincortY])
        hold on;
        %calculate overlap between amputated sides
        %overlap_indices = intersect(indices_m{3,i_task},indices_m{4,i_task});
        
        %round to nearest hundredths
        amp_pre_locations = newtargetL_m{3,i_task};
        amp_post_locations = newtargetL_m{4,i_task};
        
        %find where mean strengths are above threshold
        
        %only plots active sources with mean higher than threshold in 100ms window
%         threshold_pre = threshold_m{3,i_task};
%         threshold_post = threshold_m{4,i_task};
%         W_pre = W_m{3,i_task};
%         W_post = W_m{4,i_task};
%         active_mean_sources_pre = W_pre(W_pre>threshold_pre);
%         active_mean_sources_post = W_post(W_post>threshold_post);
%         
%         active_amp_pre_locations = amp_pre_locations(W_pre>threshold_pre,:);
%         active_amp_post_locations = amp_post_locations(W_post>threshold_post,:);
        
        %plots all active sources within 100ms window
        active_amp_pre_locations = amp_pre_locations(indices_m{3,i_task},:);
        active_amp_post_locations = amp_post_locations(indices_m{4,i_task},:);
        
        round_amp_pre_locations = round(100*active_amp_pre_locations)/100;
        round_amp_post_locations = round(100*active_amp_post_locations)/100;
        
        [overlap_indices,ia,ib] = intersect(round_amp_pre_locations(:,1),round_amp_post_locations(:,1));
        
        unique_pre_loc_ind = setdiff(1:length(active_amp_pre_locations(:,1)),ia);
        unique_post_loc_ind = setdiff(1:length(active_amp_post_locations(:,1)),ib);
        
        %plot amputated side pre and post
        plot3(active_amp_pre_locations(unique_pre_loc_ind,1),active_amp_pre_locations(unique_pre_loc_ind,2),active_amp_pre_locations(unique_pre_loc_ind,3),'Color','r','Markersize',12,'LineStyle','.');
        plot3(active_amp_post_locations(unique_post_loc_ind,1),active_amp_post_locations(unique_post_loc_ind,2),active_amp_post_locations(unique_post_loc_ind,3),'Color','b','Markersize',12,'LineStyle','.');
        
        %plot intact side
        active_int_locations = newtargetL_m{2,i_task};
        active_int_locations = active_int_locations(indices_m{2,i_task},:);
        plot3(active_int_locations(:,1),active_int_locations(:,2),active_int_locations(:,3),'Color','c','Markersize',12,'LineStyle','.');
        
        %plot overlap between amputated sides
        plot3(active_amp_post_locations(ib,1),active_amp_post_locations(ib,2),active_amp_post_locations(ib,3),'Color','g','Markersize',12,'LineStyle','.');
        %fix view
        view(0,90)
        axis on

        %create grid for show purposes
        %grid on
        axis xy

        axis equal
        axis([cdr.mincortX cdr.maxcortX cdr.mincortY cdr.maxcortY])
        axis tight
    end
            
    %plot COGs with ROI
    figure(10)
    set(gcf,'Color',[1 1 1]);
    for i_task = 1:length(grandlist)
        
        %subplot(1,length(grandlist),i_task)
        subplot(1,4,i_task)
        
        imshow(cdr.imA_resize2,'XData',[cdr.mincortX cdr.maxcortX],'YData',[cdr.maxcortY cdr.mincortY])
        hold on;
        
        %show single COGs
                
        %plot single COG
        single_mean_loc_pre = single_mean_loc_m{3,i_task};
        single_mean_loc_post = single_mean_loc_m{4,i_task};
        mirror_intact_mean_loc = mirror_single_mean_loc_m{1,i_task};
        
        plot(single_mean_loc_pre(1),single_mean_loc_pre(2),'Markersize',6,'color','r','MarkerFaceColor','r','LineStyle','s','Linewidth',3)
        plot(single_mean_loc_post(1),single_mean_loc_post(2),'Markersize',24,'color','b','MarkerFaceColor','b','LineStyle','.','Linewidth',3)
        
        plot(mirror_intact_mean_loc(1),mirror_intact_mean_loc(2),'Markersize',6,'color','c','MarkerFaceColor','c','LineStyle','^','Linewidth',3)
        
        %plot region of interest
        rot_new_perim = cdr.rot_new_perim;
        %remove points either left or right of center line
        if isequal(sides{2},'right')
            %remove right side
            rot_new_perim = rot_new_perim(find(rot_new_perim(:,1)<=0),:);
        else
            %remove left side
            rot_new_perim = rot_new_perim(find(rot_new_perim(:,1)>=0),:);
        end
        
        %draw midline between 2 points closest to edge
        [y_sort,i_sort] = sort(abs(rot_new_perim),1);
        
        jump_ycoords = diff(rot_new_perim(i_sort(:,1),2));
        jump_ind_first = find(jump_ycoords>5);
        jump_ind_first = jump_ind_first(1)+1;
        
        midline_coord = rot_new_perim(i_sort(1,1),2):rot_new_perim(i_sort(jump_ind_first,1),2);
        midline_coord = [zeros(length(midline_coord),1) midline_coord'];
        
        plot([rot_new_perim(:,1);midline_coord(:,1)],[rot_new_perim(:,2);midline_coord(:,2)],'LineWidth',7,'color','k','Linestyle','.');
        
        %fix view
        view(0,90)
        axis on

        %create grid for show purposes
        %grid on
        axis xy

        axis equal
        axis([cdr.mincortX cdr.maxcortX cdr.mincortY cdr.maxcortY])
        axis tight
    end
    
    %plot strengths and sizes of activity with respect to amount of EMG activation for each task
    figure(13)
    set(gcf,'Color',[1 1 1]);
    for i_task = 1:length(grandlist)
        subplot(1,2,1)
        hold on;
        
        %get amount of EMG activation and plot vs strengths and size
        plot(totalEMG_m(2,i_task),mean_numsources(2,i_task),'Markersize',6,'color','c','MarkerFaceColor','c','LineStyle','^','Linewidth',3);
        plot(totalEMG_m(3,i_task),mean_numsources(3,i_task),'Markersize',6,'color','r','MarkerFaceColor','r','LineStyle','s','Linewidth',3);
        plot(totalEMG_m(4,i_task),mean_numsources(4,i_task),'Markersize',24,'color','b','MarkerFaceColor','b','LineStyle','.','Linewidth',3);
        
        subplot(1,2,2)
        hold on;
        
        %get amount of EMG activation and plot vs strengths and size
        plot(totalEMG_m(2,i_task),mean_meanstrengths(2,i_task),'Markersize',6,'color','c','MarkerFaceColor','c','LineStyle','^','Linewidth',3);
        plot(totalEMG_m(3,i_task),mean_meanstrengths(3,i_task),'Markersize',6,'color','r','MarkerFaceColor','r','LineStyle','s','Linewidth',3);
        plot(totalEMG_m(4,i_task),mean_meanstrengths(4,i_task),'Markersize',24,'color','b','MarkerFaceColor','b','LineStyle','.','Linewidth',3);
    end
    %draw regression line and display regression statistics
    X_emg = totalEMG_m(2:4,:);
    X_emg = reshape(X_emg,size(X_emg,1)*size(X_emg,2),1);
    Y_numsources = mean_numsources(2:4,:);
    Y_numsources = reshape(Y_numsources,size(Y_numsources,1)*size(Y_numsources,2),1);
    Y_meanstrengths = mean_meanstrengths(2:4,:);
    Y_meanstrengths = reshape(Y_meanstrengths,size(Y_meanstrengths,1)*size(Y_meanstrengths,2),1);
    
    %regression
    if strfind(grandlist{i_task},'AK')
        P_numsources = polyfit(X_emg(2:end),Y_numsources(2:end),1);
        P_meanstrengths = polyfit(X_emg(2:end),Y_meanstrengths(2:end),1);
        
        xs = linspace(min(X_emg(2:end)),max(X_emg(2:end)),50);
    else
        P_numsources = polyfit(X_emg,Y_numsources,1);
        P_meanstrengths = polyfit(X_emg,Y_meanstrengths,1);
        
        xs = linspace(min(X_emg),max(X_emg),50);
    end
    
    ys_numsources = P_numsources(1)*xs + P_numsources(2);
    ys_meanstrengths = P_meanstrengths(1)*xs + P_meanstrengths(2);
    
    if strfind(grandlist{i_task},'AK')
        [b,bint,r,rint,stats] = regress(Y_numsources(2:end),[X_emg(2:end) ones(length(X_emg(2:end)),1)]);
        [b2,bint2,r2,rint2,stats2] = regress(Y_meanstrengths(2:end),[X_emg(2:end) ones(length(X_emg(2:end)),1)]);
    else
        [b,bint,r,rint,stats] = regress(Y_numsources,[X_emg ones(length(X_emg),1)]);
        [b2,bint2,r2,rint2,stats2] = regress(Y_meanstrengths,[X_emg ones(length(X_emg),1)]);
    end

    subplot(1,2,1)
    plot(xs,ys_numsources,'k','Linewidth',2);
    xlabel('sum of EMGs')
    ylabel('size')
    r2_numsources = stats(1)
    text(xs(end),ys_numsources(end-5),['r^2 = ',num2str(round(100*r2_numsources)/100)]);
    axis square
    subplot(1,2,2)
    plot(xs,ys_meanstrengths,'k','Linewidth',2);
    xlabel('sum of EMGs')
    ylabel('mean strengths')
    r2_meanstrengths = stats2(1)
    text(xs(end),ys_meanstrengths(end-5),['r^2 = ',num2str(round(100*r2_meanstrengths)/100)]);
    axis square
    


    dmean_total = mean(dmean_m2,1);
    numzeros = sum(dmean_m2(:,1)==0);
    dmean_total = dmean_total*length(dmean_m2(:,1))/(length(dmean_m2(:,1))-numzeros)
    
    dstd_total = std(dmean_m2,1);
    dstd_total = dstd_total*sqrt(length(dmean_m2(:,1)))/sqrt(length(dmean_m2(:,1))-numzeros)
    
    xshiftmean_total = mean(xshiftmean_m2,1);
    xshiftstd_total = std(xshiftmean_m2,1);
    yshiftmean_total = mean(yshiftmean_m2,1);
    yshiftstd_total = std(yshiftmean_m2,1);

    dmean_single_total = mean(dmean_single_m2,1);
    dstd_single_total = std(dmean_single_m2,1);
    xshiftmean_single_total = mean(xshiftmean_single_m2,1);
    xshiftstd_single_total = std(xshiftmean_single_m2,1);
    yshiftmean_single_total = mean(yshiftmean_single_m2,1);
    yshiftstd_single_total = std(yshiftmean_single_m2,1);

    results(i_k).dmean_m = dmean_m;
    results(i_k).dmean_m2 = dmean_m2;
    results(i_k).dmean_total = dmean_total';
    results(i_k).dstd_total = dstd_total';
    results(i_k).xshiftmean_m = xshiftmean_m;
    results(i_k).xshiftmean_m2 = xshiftmean_m2;
    results(i_k).xshiftmean_total = xshiftmean_total;
    results(i_k).xshiftstd_total = xshiftstd_total;
    results(i_k).yshiftmean_m = yshiftmean_m;
    results(i_k).yshiftmean_m2 = yshiftmean_m2;
    results(i_k).yshiftmean_total = yshiftmean_total;
    results(i_k).yshiftstd_total = yshiftstd_total;

    results(i_k).dmean_single_m = dmean_single_m;
    results(i_k).dmean_single_m2 = dmean_single_m2;
    results(i_k).xshiftmean_single_m = xshiftmean_single_m;
    results(i_k).xshiftmean_single_m2 = xshiftmean_single_m2;
    results(i_k).yshiftmean_single_m = yshiftmean_single_m;
    results(i_k).yshiftmean_single_m2 = yshiftmean_single_m2;
    results(i_k).dmean_single_total = dmean_single_total;
    results(i_k).dstd_single_total = dstd_single_total;
    results(i_k).xshiftmean_single_total = xshiftmean_single_total;
    results(i_k).xshiftstd_single_total = xshiftstd_single_total;
    results(i_k).yshiftmean_single_total = yshiftmean_single_total;
    results(i_k).yshiftstd_single_total = yshiftstd_single_total;

    results(i_k).single_mean_loc_m = single_mean_loc_m;
    results(i_k).mirror_single_mean_loc_m = mirror_single_mean_loc_m;
    results(i_k).ctrs_matrix = ctrs_matrix;
    results(i_k).indices_m = indices_m;
    results(i_k).mean_numsources = mean_numsources;
    results(i_k).mean_meanstrengths = mean_meanstrengths;
    results(i_k).dmean_single_amp = dmean_single_amp;
    %----------------------------------------------
    %perform NORMALIZATION
    fprintf('paused for next value of k (press enter to continue)...\n');
    %pause
end















%--------------------------------------------------------------------------
%SUB-FUNCTIONS

function [d, xshift, yshift, final_ind_m] = distlistfun(ctrs1, ctrs2, limitflag)

if limitflag == 0
    threshold = 1000;
else
    threshold = limitflag;
end

[n,p] = size(ctrs1);

d = [];
xshift = [];
yshift = [];
final_ind_m = [];
for i=1:n
    d1 = sqrt(sum((ctrs1-repmat(ctrs2(i,:),n,1)).^2,2));
    xshift1 = repmat(ctrs2(i,1),n,1) - ctrs1(:,1);     %subtract first group of centers from second group
    yshift1 = repmat(ctrs2(i,2),n,1) - ctrs1(:,2);     %subtract first group of centers from second group
    
    d2 = d1(d1<threshold);
    
    ind = (d1<threshold).*(1:n)';
    ind = ind(ind>0);
    
    if ~isempty(ind)
        ind_m = [ind repmat(i,length(ind),1)];

        d = [d; d2];
        xshift = [xshift; xshift1(ind)];
        yshift = [yshift; yshift1(ind)];
        final_ind_m = [final_ind_m; ind_m];
    end
end






