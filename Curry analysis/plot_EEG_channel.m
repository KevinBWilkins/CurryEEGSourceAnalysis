function plot_EEG_channel (subjects, channelNums)

taskLists={'abd';'ef'};
Base_dir=['F:\data\export\'];
t = (0:640-1)/256;
t = t-2;
figure (5)
hold on
for i_task = 1: length(taskLists)
    cur_task = taskLists{i_task};
    for k=1:length(subjects)
        subjectName=subjects{k};
        eeg_file=[Base_dir,subjectName,'\',cur_task,'.dat'];
        
        EEG= load (eeg_file);
        EEG_channel(k,:,:) = EEG(channelNums,:); %subject, task, channel, time
    end
    mean_EEG_channel = -1*squeeze(mean (EEG_channel,1));
    ste_EEG_channel = std (EEG_channel,1)/sqrt (length(subjects)-1);
    
    subplot(2,1,i_task)
    hold on
    errorbar(t,mean_EEG_channel,ste_EEG_channel, 'c:','LineWidth', 0.1)
    plot (t,mean_EEG_channel,'b','LineWidth', 2)
    xlabel ('Time (s)')
    ylabel ('Votage (uV)')
end

