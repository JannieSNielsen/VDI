function peak=extractPeaks_data(Xdata,Ydata)

%%extractPeaks

plot_1=0;
thres   = 1.4; % thres     = Threshold value - number of std above the mean value.
septime =  10; % septime   = Minimum Seperation time between peaks.

% PEAK OVER THRESHOLD ANALYSIS
resp = Ydata;
time = Xdata;
[peak,timepeak,error] = POT(resp,time,thres,septime);
if error == 0
   % PLOT OF RESPONSE AND PEAKS
    if plot_1==1
        figure
        plot(time,resp,'-b')
        hold on
        set(gca,'fontname','Times New Roman')
        plot(timepeak,peak,'or')
        plot([min(time) max(time)],[mean(Ydata)+thres*std(Ydata), mean(Ydata)+thres*std(Ydata)],'-k')
        xlabel('Time [sec.]','fontname','Times New Roman');
        ylabel('Response [kNm]','fontname','Times New Roman');
    end
end



