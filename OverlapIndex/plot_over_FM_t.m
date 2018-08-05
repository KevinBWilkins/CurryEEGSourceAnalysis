function plot_over_FM_t (FM, overlap_t)
time_c=[-0.6727	-0.622	-0.5713	-0.5206	-0.4699	-0.4192	-0.3685	-0.3178	-0.2671	-0.2164	-0.1657	-0.115];

y=FM;
Time_sample = size (overlap_t,2);
figure
hold on
for t=6:Time_sample
    x=overlap_t(:,t);
    p=polyfit(x,y,1);
    f=polyval(p,x);
    subplot(2,3,t-5);
    plot(x, y, 'o',x,f,'-');
    axis([0.09 0.2 0 65])
    xlabel ('Overlap')
    ylabel ('FMS')
    title (['Time =',num2str(time_c(t)),' (ms); k=', num2str(p(1))])
end
