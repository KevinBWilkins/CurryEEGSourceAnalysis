function P2Pset_distance = P2Pset_dis (point, point_set)
tmp3 = point_set-repmat(point,size(point_set,1),1);
P2Pset_distance = sqrt(tmp3(:,1).^2+tmp3(:,2).^2+tmp3(:,3).^2);
