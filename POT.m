function [peak,timepeak,error] = POT(resp,time,thres,septime)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS FUNCTION PERFORMS A PEAK OVER THRESHOLD (POT) ANALYSIS FOR A
% DATASERIES. THE RESULT IS GIVEN AS PEAKS AND THE TIME THEY OCCUR.
%
% INPUT PARAMETER:
% resp     = Response for which the POT-analysis is performed, the values
%            are given as a coulumn vector.
% time     = timeseries for response, the values are given as a coulumn
%            vector.
% thres    = Threshold value - number of std above the mean value.
% septime  = Minimum Seperation time between peaks.
%
%
% OUTPUT PARAMETERS:
% peak     = Coulumn vector which contains the extrated peaks
% timepeak = Coulumn vector which contains the time for the individual
%            peaks
% error    = 0 = At lest one peak in response
%            1 = No Peaks in response
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATION OF THRESHOLD VALUE
my  = mean(resp);
sig = std(resp);
threshold = my + thres * sig;

% CALCULATION OF LENGTH OF RESPONSE VECTOR AND STARTING COUNTER FOR NUMBER
% OF PEAKS K
n = length(resp);
k=0;

% CASE 1: IF THE FIRST VALUE IS ABOVE THRESHOLD !
if(resp(1)>threshold)
    for j=1:n-1
        if resp(j) > threshold & resp(j+1) < threshold
            up   = 1;
            down = j+1;
            break
        end
        up   = j;
        down = n;
    end
    k=k+1;

    [peaks(k), index] = max(resp(up:down));
    timepeaks(k) = time(up-1+index);
end

% CASE 2: PEAKS FOR REST OF RESPONSE VECTOR OR FIRST VALUE BELOW THRESHOLD
for i=1:n-1
    if resp(i) < threshold & resp(i+1) > threshold
        for j=i+1:n-1
            if resp(j) > threshold & resp(j+1) < threshold
                up   = i;
                down = j+1;
                break
            end
            up   = i;
            down = n;
        end
        k=k+1;

        [peaks(k), index] = max(resp(up:down));
        timepeaks(k) = time(up-1+index);
    end
end

% SPLITTING OF THE PEAKS SO THEY ARE INDEPENDENT - SEPERATION TIME
% "septime"
if k > 0
    i=0;
    for inum=1:n
        i=i+1;

        space=0;
        position=0;
        position=zeros(1,k-1);
        for h=1:k-1
            space(h)=timepeaks(h+1)-timepeaks(h);
            if space(h) < septime
                position(h)=1;
            end
        end

        if position(i) == 1
            low=i;
            for j=i:k-1
                if position(j) == 0
                    high=j;
                    break
                end
                high=k-1;
            end

            [value, location] = max(peaks(low:high));

            %       MAX FIRST IN PEAK SERIES - DELETE VALUE AFTER PEAK
            if location == 1
                peaks(low+1:k-1)=peaks(low+2:k);
                timepeaks(low+1:k-1)=timepeaks(low+2:k);
                k=k-1;
                i=i-1;
            end

            %       MAX MIDDLE IN PEAK SERIES - DELETE VALUE BEFORE AND AFTER PEAK
            if location > 1 & location < high-low+1
                peaks(low+location-2)=peaks(low+location-1);
                peaks(low+location-1:k-2)=peaks(low+location+1:k);
                timepeaks(low+location-2)=timepeaks(low+location-1);
                timepeaks(low+location-1:k-2)=timepeaks(low+location+1:k);
                k=k-2;
                i=i-1;
            end

            %       MAX END IN PEAK SERIES - DELETE VALUE BEFORE PEAK
            if location == high-low+1 & location > 1
                peaks(high-1:k-1)=peaks(high:k);
                timepeaks(high-1:k-1)=timepeaks(high:k);
                k=k-1;
                i=i-1;
            end
        end
        if i == k-1
            break
        end
    end
    peak=peaks(1:k);
    timepeak=timepeaks(1:k);
    error=0;
elseif k == 0
    peak     = 0;
    timepeak = 0;
    error    = 1;
end