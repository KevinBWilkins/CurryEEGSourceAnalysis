function plot_act_time_loc_S1 (act_ratio_loc_time, LI)

[K, I, loc_no, time_no]= size(act_ratio_loc_time);

winOverlap=0;
winSize=1;
timePoints=65;
POINT_T_DELTA        =  0.49;
start_T = 10;
taskList ={'m'};
k=1;

X=0:POINT_T_DELTA:(timePoints-1)*POINT_T_DELTA;
X=X+start_T;
X=fliplr(X);
%     Y=[-10:-1,1:10];
%     Y=1:18;

    

for i=1:I
    figure (5+i); clf;        
    subplot(2,3,2);
    %         X=size(act_ratio_loc_time,4);
    nY=size(act_ratio_loc_time,3); Y=1:nY;
    Y = [-1*nY/2:-1,1:nY/2];
    %         surf(X,Y,squeeze(act_ratio_loc_time(k,i,:,:)));
    tmp=squeeze(act_ratio_loc_time(k,i,:,:));
%     if i>0
%         tmp=flipud(tmp);
%     end
    
    imagesc(Y,X,flipud(flipud(tmp)')); %colorbar;

    %         set(gca,'YTickLabel',mat2cell ([10:40]))
    %         set(gca,'YTickLabel',{'lip','lip-hand-R','hand-R','elbow-R','shoulder-R','shoulder-L','elbow-L','hand-L','lip-hand-L','lip'})
    % title (taskList{i})
    %         view (-62,83);


    tim = sum (tmp, 1);
    loc = sum (tmp, 2);
    loc=flipud(loc);

    [b,a] = butter(9,500*POINT_T_DELTA*0.001,'low');
    tim_evn = filtfilt(b,a,tim);
   
    
    
    subplot(2,3,1)
    hold on
    plot(tim); plot(tim_evn, 'g');
    view (-90, -90)

    subplot(2,3,3)
    hold on
    plot(squeeze(LI(k,i,:)));
    plot(zeros(1,timePoints), 'r:');
    view (90, 90)

    subplot(2,3,5)
    plot(loc)
    view (0, -90)
end
        
   
