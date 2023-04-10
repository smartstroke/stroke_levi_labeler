clear board; clear all;
% % % Initialize the Arduino board % % %

board = serialport("COM5",115200);
configureTerminator(board,"CR/LF");
flush(board);

% % % Received data variables % % %
data0 = 0;
data1 = 0;
time = 0;

% % % Scalars % % %
peak_tot_lag = 0; % Total time difference calculated from location of peaks
trough_tot_lag = 0; % Total time difference calculated from location of troughs
dt_peak_tot_lag = 0; % Derivitive of pressure waveform Total time lag calculated from peaksfind_curveBase1_bool_min
dt_trough_tot_lag = 0; % Derivitive of pressure waveform Total time lag calculated from troughs
sync = 0; % Variable that is the average of all 'lag' variables
tol = 0.05; % Tolerance to scale max peak height;
curr_num_peaks0 = 0; % Compares with prev_data0,1_peaks
curr_num_peaks1 = 0;

% % % Loop control % % %
duration = 10; % duration of while loop in sec
stepTime = 0.001; % Used for pause()

% % % Comparitors % % % Holds number of peaks, and peaks of derivative
prev_data0_peaks = 0; 
prev_data1_peaks = 0;
prev_data0_dpeaks = 0;
prev_data1_dpeaks = 0;

% % % Indices % % %
dataIndex = 1;
locIndex = 1; % Stores Peak location Index
dlocIndex = 1; % Stores Derivative Peak location index
indexHolder = 0;

% % % Array Definitions % % %
peaks_lag = zeros(1,2000*duration); % Array to track lag detemined from peak location at each index (peaks_lag)
toughs_lag = zeros(1,2000); % Array to track lag detemined from trough location at each index (troughs_lag)
dt_peaks_lag = zeros(1,2000); % Array to track lag detemined from peak location at each index (dt_peaks_lag)
find_curveBase_bool = zeros(1,2000*duration); % Logical arrays
find_curveBase1_bool = zeros(1,2000*duration);
% Store data from sensors
A0Array = zeros(1,2000*duration);
A1Array = zeros(1,2000*duration);
Tstamp = zeros(1,2000*duration);

% % % Initialize live plot % % %

subplot(2, 1, 1);
hold on
line0 = animatedline('Color','r','MaximumNumPoints',1000000);
line1 = animatedline('Color','b','MaximumNumPoints',1000000);
% xlim([0,duration]);
title("Recorded Data");
xlabel("Time (millisec)");
ylim([0,5000]);
legend('Signal 0', 'Signal 1');

subplot(2, 1, 2);
line2 = animatedline('Color','green','MaximumNumPoints',1000000);% Used to be time_diff
title("Time Lag");
xlabel("Time");
% subplot(3, 1, 3)
% line3 = line(start-start,dt_peak_tot_lag,'Color','magenta');
% title("Time Lag from Derivative");

flush(board);
pause(2);
while (time<=duration*1000)
    board
    % Don't allow buffer to stack up
%     if board.NumBytesAvailable>100
%         flush(board);
%     end
%     write(board,'a','char');
    if board.NumBytesAvailable>120
        indexHolder = dataIndex;
        while(board.NumBytesAvailable>120)
            data = read(board,3,'uint32'); 
            data0 = data(1);
            data1 = data(2);
            time = data(3);
            A0Array(indexHolder) = data0;
            Tstamp(indexHolder) = time/1000;
            A1Array(indexHolder) = data1;
            indexHolder = indexHolder+1;
        end
    else
        %Sending three 4-byte long ints so read as 3 uint32
        data = read(board,3,'uint32'); 
        data0 = data(1);
        data1 = data(2);
        time = data(3);
    
        A0Array(dataIndex) = data0;
        Tstamp(dataIndex) = time/1000;
        A1Array(dataIndex) = data1;

    end
    

    

    deriv0 = gradient(A0Array);
    deriv1 = gradient(A1Array);

    [dpeaks0, dploc0] = findpeaks(deriv0,"MinPeakDistance",60,'MinPeakHeight',max(deriv0)-3*std(deriv0));%,'MinPeakProminence',1.5);
    [dpeaks1, dploc1] = findpeaks(deriv1,"MinPeakDistance",60,'MinPeakHeight',max(deriv1)-3*std(deriv1));%,'MinPeakProminence',1.5);
        
    % maybe get ndpeaks from the width of peak + dploc?
    % ex: ndploc0 = dploc0 + Tstamp(find(Tstamp == w0));
%     [ndpeaks0, ndploc0] = findpeaks(-deriv0,"MinPeakDistance",60);%,'MinPeakHeight',max(-deriv0)-3*std(-deriv0),'MinPeakProminence',1);%,'MinPeakProminence',0.1);
%     [ndpeaks1, ndploc1] = findpeaks(-deriv1,"MinPeakDistance",60);%,'MinPeakHeight',max(-deriv1)-3*std(-deriv1),'MinPeakProminence',1);%,'MinPeakProminence',0.1);



%     % Find base edges of curves
    find_curveBase0_bool_max = (tol*max(A0Array)<A0Array);
    find_curveBase0_bool_min = (A0Array>min(A0Array));
    find_curveBase0_bool = find_curveBase0_bool_max.*find_curveBase0_bool_min;
    find_curveBase0_bool = [false,find_curveBase0_bool,false];
    start0 = strfind(find_curveBase0_bool,[0 1]);
    end0 = strfind(find_curveBase0_bool,[1 0]);
    find_curveBase1_bool_max = (tol*max(A1Array)<A1Array);
    find_curveBase1_bool_min = (A1Array>min(A1Array));
    find_curveBase1_bool = find_curveBase1_bool_max.*find_curveBase1_bool_min;
    find_curveBase1_bool = [false,find_curveBase1_bool,false];
    start1 = strfind(find_curveBase1_bool,[0 1]);
    end1 = strfind(find_curveBase1_bool,[1 0]);
    start0 = start0 - 1; % -1 corrects for current loop
    start1 = start1 - 1;
    

%      %Find Peaks
    [peaks0, ploc0] = findpeaks(A0Array);%,"MinPeakDistance",6,'MinPeakHeight',max(A0Array)-std(A0Array));%,'MinPeakProminence',3);%;%max(A0Array)-2*std(A0Array),'WidthReference','halfheight');
    [peaks1, ploc1] = findpeaks(A1Array);%,"MinPeakDistance",6,'MinPeakHeight',max(A1Array)-std(A0Array));%,'MinPeakProminence',3);%2*std(A0Array),'WidthReference','halfheight');

    % Find Troughs
%     [npeaks0, nploc0] = findpeaks(-A0Array,"MinPeakDistance",4,'MinPeakHeight',max(-A0Array)-2*std(-A0Array),'MinPeakProminence',0.1);
%     [npeaks1, nploc1] = findpeaks(-A1Array,"MinPeakDistance",4,'MinPeakHeight',max(-A1Array)-2*std(-A1Array),'MinPeakProminence',0.1);
    
    % Get number of peaks
    curr_num_peaks0 = length(ploc0);
    curr_num_peaks1 = length(ploc1);
    curr_num_dpeaks0 = length(dploc0);
    curr_num_dpeaks1 = length(dploc1);

    % Skip if no peaks, Finds peak_tot_lag and sets peaks_lag(index)
    if ~isempty(ploc0)&&~isempty(ploc1)
        % Check if there are new peaks to find lag for
        if (curr_num_peaks1>prev_data1_peaks)||(curr_num_peaks0>prev_data0_peaks)
            % Check if peak loaction arrays are same size
            if length(ploc0)==length(ploc1)
                peak_tot_lag = peak_tot_lag + (Tstamp(ploc1(locIndex)) - Tstamp(ploc0(locIndex)));
                peaks_lag(dataIndex) = Tstamp(ploc1(locIndex))-Tstamp(ploc0(locIndex));
                locIndex = locIndex + 1;
            else % Take the last index for both
                peak_tot_lag = peak_tot_lag + (Tstamp(ploc1(end)) - Tstamp(ploc0(end)));
                locIndex = max(length(ploc1),length(ploc0));

            end
            prev_data0_peaks = curr_num_peaks0;
            prev_data1_peaks = curr_num_peaks1;
            % Do the same for the Troughs
%             if ~(isempty(nploc1)&&isempty(nploc0))
%                 ntlag = ntlag + (Tstamp1(nploc1(end)) - Tstamp(nploc0(end)));
%                 troughs_lag(dataIndex) = Tstamp1(nploc1(end))-Tstamp(nploc0(end));
%             end
        end
    end

%     Finds the time lag calculated from the peaks of the derivative
%     if ~isempty(dploc0)&&~isempty(dploc1)
%         if (curr_num_dpeaks1>prev_data1_dpeaks)||(curr_num_dpeaks0>prev_data0_dpeaks)
%             if (length(dploc0)||length(dploc1))==dlocIndex
% %                 peak_tot_lag = peak_tot_lag + (ploc1(end) - ploc0(end) + (dploc1(end) - dploc0(end)))/2;
%                 dt_peak_tot_lag = dt_peak_tot_lag + (Tstamp(dploc1(dlocIndex)) - Tstamp(dploc0(dlocIndex)));
%                 dt_peaks_lag(dataIndex) = Tstamp(dploc1(dlocIndex)) - Tstamp(dploc0(dlocIndex));
%                 dlocIndex = dlocIndex + 1;
%             else
%                 dt_peak_tot_lag = dt_peak_tot_lag + (Tstamp(dploc1(end)) - Tstamp(dploc0(end)));
%                 dt_peaks_lag(dataIndex) = Tstamp(dploc1(end)) - Tstamp(dploc0(end));
%                 %dt_trough_tot_lag = dt_trough_tot_lag + (Tstamp1(ndploc1(end)) - Tstamp(ndploc0(end)));
%                 dlocIndex = max(length(dploc1),length(dploc0));
%             prev_data0_dpeaks = curr_num_dpeaks0;
%             prev_data1_dpeaks = curr_num_dpeaks1;
% 
%             end
%         end
%     end
    if(indexHolder>dataIndex)
        while(indexHolder>dataIndex)
            addpoints(line0,Tstamp(dataIndex),A0Array(dataIndex));
            addpoints(line1,Tstamp(dataIndex),A1Array(dataIndex));
            addpoints(line2,Tstamp(dataIndex),peaks_lag(dataIndex));
            dataIndex = dataIndex+1;
        end
        
    else
        addpoints(line0,Tstamp(dataIndex),data0);
        addpoints(line1,Tstamp(dataIndex),data1);
        addpoints(line2,Tstamp(dataIndex),peaks_lag(dataIndex));%peak_tot_lag);

    end
    drawnow limitrate %nocallbacks
    dataIndex = dataIndex + 1;
    
end

dataIndex
peaks_lag = nonzeros(peaks_lag);
avg_of_peaks_arr = mean(peaks_lag);
max_lag = max(peaks_lag); % Maximum instantaneous lag in seconds
min_lag = min(peaks_lag);
sum_peaks_lag = sum(peaks_lag); % total Time lag in seconds should == peaks_tot_lag

% freq0 = length(peaks0)/timer % strokes (peaks) per second
% freq1 = length(peaks1)/timer

% ppos0 = length(ploc0); % total number of peaks
% ppos1 = length(ploc1);


% Show peaks on plot

% figure;
% plot(A0Array,"Color",'blue');
% hold on
% plot(Tstamp,deriv0,"Color",'red');
% plot(A1Array,"Color",'cyan')
% plot(Tstamp,deriv1,"Color",'blue')
% % plot(ploc0,peaks0,'r^')
% % plot(ploc1,peaks1,'g^')
% plot(dploc0,dpeaks0,'g^')
% plot(dploc1,dpeaks1,'r^')
% plot(ndploc0,-ndpeaks0,'g^')
% plot(ndploc1,-ndpeaks1,'r^')


% % Plot of curves with profile markings
% figure; 
% plot(Tstamp,A0Array,"Color",'red'); hold on; plot(Tstamp(ploc0),peaks0,'go');
% plot(Tstamp,A1Array,"Color",'blue'); plot(Tstamp(ploc1),peaks1,'mo');
% plot(Tstamp(dploc0),dpeaks0,'g^'); plot(Tstamp(dploc1),dpeaks1,'m^'); 
% % plot(Tstamp(ndploc0),ndpeaks0,'gv'); plot(Tstamp(ndploc1),ndpeaks1,'mv'); 
% plot(Tstamp(pstart0),A0Array(pstart0),'g>',Tstamp(pend0),A0Array(pend0),'g<');
% plot(Tstamp(pstart1),A1Array(pstart1),'m>',Tstamp(pend1),A1Array(pend1),'m<');

% Find the sync factor, should add a length check
% for k = 1:length(ploc0)
%     sync(k) = (Tstamp1(ploc1(k)) - Tstamp(ploc0(k))) + (Tstamp1(dploc1(k)) - Tstamp(dploc0(k))) + (Tstamp1(ndploc1(k)) - Tstamp(ndploc0(k)));
% end

% prints

% avg_of_peaks_arr
% max_lag
% min_lag
% sum_peaks_lag
% peak_tot_lag
% dt_peak_tot_lag
% % dt_trough_tot_lag
% idx = find(peaks_lag==min_lag);
% time_of_min = (idx);
% time_of_min
% sum(sync)/3


clear all;
%SADC(abs(SADC-mean(SADC))<2*std(SADC)) = 0; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!