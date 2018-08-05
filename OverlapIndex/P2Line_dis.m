function allP2Line_distance = P2Line_dis(allPoints, P14line, P24line)
%     for i=1:size(allPoints,1)
%         cur_point = allPoints(i,:);
%         allP2Line_distance (i) = norm( cross ((P24line-P14line),(P14line-cur_point)))/ norm (P24line-P14line);
%     end
tmp = repmat(P24line-P14line,size(allPoints,1),1);
tmp1 = repmat(P14line,size(allPoints,1),1)-allPoints;
tmp3=cross(tmp,tmp1,2);
tmp4 = sqrt(tmp3(:,1).^2+tmp3(:,2).^2+tmp3(:,3).^2);
allP2Line_distance= tmp4./norm (P24line-P14line);