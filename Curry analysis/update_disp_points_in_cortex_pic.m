function update_disp_points_in_cortex_pic(rot_cortexL, ROI_ind, ROI_ind_group, mincortX, maxcortX, maxcortY, mincortY, imA_resize2)

ROI_cortex=rot_cortexL(ROI_ind,:);
left_ROI_cortex = rot_cortexL(ROI_ind(ROI_ind_group>0),:);
right_ROI_cortex = rot_cortexL(ROI_ind(ROI_ind_group<1),:);

figure(5)
clf

set(gcf,'Color',[1 1 1]);
imshow([mincortX maxcortX],[maxcortY mincortY],imA_resize2)
axis on
axis([mincortX maxcortX mincortY maxcortY])
axis xy
hold on;

plot3(left_ROI_cortex(:,1),left_ROI_cortex(:,2),left_ROI_cortex(:,3),'g.', 'MarkerSize',6);
plot3(right_ROI_cortex(:,1),right_ROI_cortex(:,2),right_ROI_cortex(:,3),'r.', 'MarkerSize',6);

view(0,90)
axis on
axis equal
axis([mincortX maxcortX mincortY maxcortY])

