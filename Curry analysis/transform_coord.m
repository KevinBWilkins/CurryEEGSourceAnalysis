tempdir = pwd;
cd 'C:\albertcurry\mridata\LS072006_mricoord\new\0025';

V = spm_vol('sstimpat1a-0025-00001-000192-00.img');

%sample coordinates (in voxels)
% %GM
% vox_coord = [37 145 89.5; 134.7 114.6 12.7; 138.3 108 156.2];
% real_coord = [6.3 115.5 -27.7; -73.3 20.5 -59.3; 70 11.9 -63.6];

%LSpost1
vox_coord = [48.7 147.1 92.5; 146.1 109.5 18.5; 138 111.3 178.9];
real_coord = [-3.8 120.4 -36.3; -76.8 23.1 -75.9; 83.5 31.2 -69.6];

%apply transformation (from voxel space to mm)
T = V.mat;

%T is 4x4 affine transformation
T1 = T(1:3,1:3);
inv_T1 = inv(T1);

T2 = T(1:3,4);

new_real_coord = (T1*vox_coord')' + repmat(T2',3,1)
new_vox_coord = (inv_T1*(real_coord - repmat(T2',3,1))')'

%now normalized coordinates
load 'sstimpat1a-0025-00001-000192-00_sn.mat'
% 
% %Affine is not correct normalization matrix
% A1 = Affine(1:3,1:3);
% A2 = Affine(1:3,4);

Q = VG(1).mat*inv(Affine)/VF.mat;

Q1 = Q(1:3,1:3);
inv_Q1 = inv(Q1);

Q2 = Q(1:3,4);

new_norm_real_coord = (Q1*real_coord')' + repmat(Q2',3,1)   %transforms real world points to normalized real world points

% V2 = spm_vol('wsGUILLERMOCAMRI-0003-00001-000001-00.img');
% real_coord2 = [.3 72.8 -37.1; -75.4 -20.3 -37.1; 77.1 -14.6 -44];
% vox_coord2 = [39.9 93.4 7.4; 77.7 46.9 7.4; 1.5 49.7 4.0];
% 
% T3 = V2.mat;
% T4 = T3(1:3,1:3);
% inv_T4 = inv(T4);
% T5 = T3(1:3,4);
% 
% new_real_coord2 = (T4*vox_coord2')' + repmat(T5',3,1)
% new_vox_coord2 = (inv_T4*(real_coord2 - repmat(T5',3,1))')'


%load cdr file and transform coordinates to normalized mri space

cd(tempdir)
cdr_file_name = 'C:\Albert Chen\Subject\LSpost1\new\abd_un_mri.cdr';
fprintf(1, 'get cortex locations... ')
[cortexL_mri_voxel,Lcount,LNR,TM]=read_Curry_file3_TM(cdr_file_name,'LOCATION',0,0);
%cortexL_mri_voxel are VOXEL coordinates
%TM is transformation matrix from MRI coord to curry coordinates
fprintf(1, 'done\n')

max(cortexL_mri_voxel,[],1)
min(cortexL_mri_voxel,[],1)

%curry mri coordinates are screwed up in relation to spm- 
%curry x is equal to -z with some offset,
%curry y is equal to x
%curry z is equal to y

%curry voxel space to spm voxel space transformation




xoffset = 228;  %THIS IS MOST IMPORTANT TO GET RIGHT!!! HOW TO GET THIS ACCURATELY???




cortexL_mri_spm_voxel(:,1) = cortexL_mri_voxel(:,2);
cortexL_mri_spm_voxel(:,2) = cortexL_mri_voxel(:,3);
cortexL_mri_spm_voxel(:,3) = -cortexL_mri_voxel(:,1)+xoffset;



TM2 = TM(2:4,1);    %translation
TM1 = TM(2:4,2:4);  %rotation and scaling

inv_TM1 = inv(TM1);

%cortexL_mri = (inv_TM1*(cortexL - repmat(TM2',LNR(1),1))')';
cortexL_curry = (TM1*cortexL_mri_voxel')' + repmat(TM2',LNR(1),1);
%mirror and rotate 180 degrees- negate x and y coordinates
cortexL_curry = [-cortexL_curry(:,1:2) cortexL_curry(:,3)];     %this matches curry mri locations to curry real world locations


%normalize cortexL_mri locations
cortexL_mri_spm_real = (T1*cortexL_mri_spm_voxel')' + repmat(T2',LNR(1),1); %this transforms cortex locations in spm voxel space to spm real world space
max(cortexL_mri_spm_real,[],1)
min(cortexL_mri_spm_real,[],1)

cortexL_mri_spm_real_norm = (Q1*cortexL_mri_spm_real')' + repmat(Q2',LNR(1),1); %this transforms cortex locations in spm real world space to spm normalized real world space
max(cortexL_mri_spm_real_norm,[],1)
min(cortexL_mri_spm_real_norm,[],1)


figure(5)
plot3(cortexL_mri_spm_real_norm(:,1),cortexL_mri_spm_real_norm(:,2),cortexL_mri_spm_real_norm(:,3),'LineWidth',7,'color','b','LineStyle','.')


