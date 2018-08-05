%run this following lines in find_target_PC4 before set the targetL for ef.

indtmp=find(ROI_ind>1992);
ROI_ind1=ROI_ind;
ROI_ind1(indtmp)=ROI_ind(indtmp)-1;
ROI_ind2=ROI_ind;
ROI_ind=ROI_ind1;
