%This function write the mean results back to the CDR file. So curry can
%read it and display it. Example command line:
%results2cdr('F:\data\inverse_results\JW-BCI\new\HCL-T.dat', 'F:\data\inverse_results\JW-BCI\new\newHCL-T.cdr', 'F:\data\inverse_results\JW-BCI\points\Mask','F:\data\inverse_results\JW-BCI\new\newHCL-T.dat', 0);
%cortex_ind: 0-right, 1-left, 2-all
function results2cdr(targetC_file, org_CDR_file, MaskFile,outputfile,cortex_ind);

TimePoints=20;
plotCort=1;
load (MaskFile)
ROI_ind=ROI_ind(find(ROI_ind_group==cortex_ind));
%------------------------------------
%read the new strength that need to display
fprintf(1, 'read the strength....')
newtargetC=load (targetC_file);

%----------------------------------------------------------------------    
%get all cortex locations
fprintf(1, 'get cortex locations... ')
[cortexL,targetL,targetC]=find_target_PC5(org_CDR_file,MaskFile,TimePoints,plotCort,cortex_ind);
fprintf(1, 'done\n')

%---------------------------------------------------
%creat the new strength file for different locations
fprintf(1, 'creat the new strength file for different locations... ')
[Nrow,Ncol]=size(newtargetC);
[n,m]=size(cortexL);
strength=zeros(n,Ncol);
strength(ROI_ind,:)=newtargetC;
save (outputfile,'strength','-ascii');
% breakstr=['new points'];
% for k=1:Ncol
%     cur_str=strength(:,k);
%     save (outputfile,'strength(:,k)','-ascii','-append');
%     save (outputfile,'breakstr','-ascii','-append');
% end
fprintf(1, 'done\n')