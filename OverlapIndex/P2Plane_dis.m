function allP2Plane_dis = P2Plane_dis (all_points, Plane_para)

all_points = cat(2, all_points, ones (size(all_points,1),1) );
allP2Plane_dis = abs(dot (all_points, repmat(Plane_para, size(all_points,1),1),2))./norm(Plane_para(1:3));