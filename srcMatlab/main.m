clear all;

% % % Initialize the Arduino board % % %
board = serialport("COM3",19200);
configureTerminator(board,"CR/LF");

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
duration = 30; % duration of while loop in sec
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

% % % Array Definitions % % %
peaks_lag = zeros(1,1000*duration); % Array to track lag detemined from peak location at each index (peaks_lag)
toughs_lag = zeros(1,1000); % Array to track lag detemined from trough location at each index (troughs_lag)
dt_peaks_lag = zeros(1,1000); % Array to track lag detemined from peak location at each index (dt_peaks_lag)
find_curveBase_bool = zeros(1,1000*duration); % Logical arrays
find_curveBase1_bool = zeros(1,1000*duration);

% Store data from sensors
A0Array = zeros(1,1000*duration);
A1Array = zeros(1,1000*duration);
Tstamp = zeros(1,1000*duration);

% % % Initialize live plot % % %
subplot(2, 1, 1);
hold on
line0 = animatedline('Color','r','MaximumNumPoints', 1000000);
line1 = animatedline('Color','b','MaximumNumPoints', 1000000);
title("Recorded Data");
xlabel("Time (seconds)");
x1 = 0;
dx = 10;
x=[0:100000];
y1 = 0;
y2 = 5000;
axis([x1 x1+dx y1 y2]);
legend('Signal 0', 'Signal 1');

% subplot(2, 1, 2);
% line2 = animatedline('Color','green','MaximumNumPoints',1000000);% Used to be time_diff
% title("Time Lag");
% xlabel("Time");

flush(board);
pause(0.01);
ticStart = tic;
time0 = 0;
while (toc(ticStart)<= duration)
    %board
    [time, data0, data1] = fGetData(board, false);

    if (dataIndex == 1)
        time0 = time;
    end

    A0Array(dataIndex) = data0;
    Tstamp(dataIndex) = ((time + 1) - time0)/1000;
    A1Array(dataIndex) = data1;

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

    addpoints(line0,Tstamp(dataIndex),data0);
    addpoints(line1,Tstamp(dataIndex),data1);

    xlim([Tstamp(dataIndex)-dx, Tstamp(dataIndex)+(dx/2)])

%     addpoints(line2,Tstamp(dataIndex),peak_tot_lag);

    drawnow limitrate %nocallbacks

    dataIndex = dataIndex + 1;
end

dataIndex;
peaks_lag = nonzeros(peaks_lag);
avg_of_peaks_arr = mean(peaks_lag);
max_lag = max(peaks_lag); % Maximum instantaneous lag in seconds
min_lag = min(peaks_lag);
sum_peaks_lag = sum(peaks_lag); % total Time lag in seconds should == peaks_tot_lag

% print
avg_of_peaks_arr;
max_lag;
min_lag;
sum_peaks_lag;
peak_tot_lag;
dt_peak_tot_lag;
% dt_trough_tot_lag
idx = find(peaks_lag==min_lag);
time_of_min = (idx);
time_of_min;
% sum(sync)/3


clear board;
%SADC(abs(SADC-mean(SADC))<2*std(SADC)) = 0; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
