clear board; clear all;
% % % Initialize the Arduino board % % %
port_list = serialportlist;
board = serialport(port_list(end),115200);
configureTerminator(board,"CR/LF");
flush(board);

% % % Received data variables % % %
data0 = 0;
data1 = 0;
time = 0;

% peaks0 = [0];
avg_width0 = 0;
% w0 = ones(1,1);
meanCycle0 = 0;
meanCycle1 = 0;

% % % Scalars % % %
tot_lag_peaks = 0; % Total time difference calculated from location of peaks
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

% % % Array Definitions % % %
peaks_lag = zeros(1,15*duration); % Array to track lag detemined from peak location at each index (peaks_lag)
toughs_lag = zeros(1,150); % Array to track lag detemined from trough location at each index (troughs_lag)
dt_peaks_lag = zeros(1,150); % Array to track lag detemined from peak location at each index (dt_peaks_lag)
find_curveBase_bool = zeros(1,15*duration); % Logical arrays
find_curveBase1_bool = zeros(1,15*duration);
total_lag_peaks_arr = zeros(1,15*duration);
% Store data from sensors
A0Array = zeros(1,15*duration);
A1Array = zeros(1,15*duration);
Tstamp = zeros(1,15*duration);

% % % Initialize live plot % % %

plt_color = '';

p1 = subplot(2, 1, 1);
hold on
line0 = animatedline('Color','r','MaximumNumPoints',1000000);
line1 = animatedline('Color','b','MaximumNumPoints',1000000);
line0p = animatedline('Color','green','Marker','o','LineStyle','none');
line1p = animatedline('Color','m','Marker','o','LineStyle','none');
% xlim([0,duration]);
title("Live Pressure Data");
ylim([0,5000]);
legend('Signal 0', 'Signal 1');
ylabel(p1,'Amplitude (mV)')
dx = 10;

p2 = subplot(2, 1, 2);
line2 = animatedline('Color','green','MaximumNumPoints',1000000, 'LineWidth',2);% Used to be time_diff
line2r = animatedline('Color','green','Marker','o','LineWidth',2,'LineStyle','none');
ylabel(p2,'Time Difference')
xlabel("Time (sec)");
color_tracer = 0;

% subplot(3, 1, 3)
% line3 = line(start-start,dt_peak_tot_lag,'Color','magenta');
% title("Time Lag from Derivative");
linkaxes([p1 p2],'x')
loop_counter = 1;
pause(2);
flush(board);
while (true)%(time<=duration*1000)
    %board
    % Don't allow buffer to stack up
    if board.NumBytesAvailable>72
        flush(board);
    end

    write(board,'a','char');

    %Sending three 4-byte long ints so read as 3 uint32

    data = read(board,3,'uint32'); 
    
    data0 = data(1);
    data1 = data(2);
    time = data(3);
    A0Array(dataIndex) = data0;
    Tstamp(dataIndex) = time/1000;
    A1Array(dataIndex) = data1;

    deriv0 = gradient(A0Array);
    deriv1 = gradient(A1Array);

    
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
    

%      %Find Peaks floor((loop_counter/duration)+((prev_data0_peaks+1)*19) ((time/1000)/(avg_width0+1))
    [peaks0, ploc0,w0,~] = findpeaks(A0Array,"MinPeakDistance",duration,'MinPeakHeight',max(A0Array)/3);%,'MinPeakProminence',2);%;%max(A0Array)-2*std(A0Array),'WidthReference','halfheight');
    [peaks1, ploc1] = findpeaks(A1Array,"MinPeakDistance",duration,'MinPeakHeight',max(A1Array)/3);%,'MinPeakProminence',2);%2*std(A0Array),'WidthReference','halfheight');

    % Find Troughs
%     [npeaks0, nploc0] = findpeaks(-A0Array,"MinPeakDistance",4,'MinPeakHeight',max(-A0Array)-2*std(-A0Array),'MinPeakProminence',0.1);
%     [npeaks1, nploc1] = findpeaks(-A1Array,"MinPeakDistance",4,'MinPeakHeight',max(-A1Array)-2*std(-A1Array),'MinPeakProminence',0.1);
    [dpeaks0, dploc0,p0] = findpeaks(deriv0,"MinPeakDistance",duration,'MinPeakHeight',max(deriv0)/2);%,'MinPeakProminence',1.5);
    [dpeaks1, dploc1,p1] = findpeaks(deriv1,"MinPeakDistance",duration,'MinPeakHeight',max(deriv1)/2);%,'MinPeakProminence',1.5);
        
    % maybe get ndpeaks from the width of peak + dploc?
    % ex: ndploc0 = dploc0 + Tstamp(find(Tstamp == w0));
    [ndpeaks0, ndploc0] = findpeaks(-deriv0,"MinPeakDistance",duration,'MinPeakHeight',max(-deriv0)/2);%,'MinPeakHeight',max(-deriv0)-3*std(-deriv0),'MinPeakProminence',1);%,'MinPeakProminence',0.1);
    [ndpeaks1, ndploc1] = findpeaks(-deriv1,"MinPeakDistance",duration,'MinPeakHeight',max(-deriv1)/2);%,'MinPeakHeight',max(-deriv1)-3*std(-deriv1),'MinPeakProminence',1);%,'MinPeakProminence',0.1);



    % Get number of peaks
    curr_num_peaks0 = length(ploc0);
    curr_num_peaks1 = length(ploc1);
    curr_num_dpeaks0 = length(dploc0);
    curr_num_dpeaks1 = length(dploc1);
    peaks_lag(dataIndex) = 0;
%     % Skip if no peaks, Finds tot_lag_peaks and sets peaks_lag(index)
    if ~isempty(ploc0)&&~isempty(ploc1)
        disp("Level 0")
        if (curr_num_peaks1>prev_data1_peaks)||(curr_num_peaks0>prev_data0_peaks)
            disp("Level 1")
            if (length(ploc0)&&length(ploc1))==locIndex
                disp("Level 2A")
%                 tot_lag_peaks = tot_lag_peaks + (ploc1(end) - ploc0(end) + (dploc1(end) - dploc0(end)))/2;
                tot_lag_peaks = tot_lag_peaks + (Tstamp(ploc1(locIndex)) - Tstamp(ploc0(locIndex)));
%                 total_lag_peaks_arr(dataIndex) = tot_lag_peaks;
                peaks_lag(dataIndex) = Tstamp(ploc1(locIndex)) - Tstamp(ploc0(locIndex));
                locIndex = locIndex + 1;
            else%if length(ploc0)>length(ploc(1)
                if abs(Tstamp(ploc0(end))-Tstamp(ploc1(end)))>3
                    peaks_lag(dataIndex) = 0;%peaks_lag(dataIndex-1);
%                     total_lag_peaks_arr(dataIndex) = total_lag_peaks_arr(dataIndex-1);
                else
                    tot_lag_peaks = tot_lag_peaks + (Tstamp(ploc1(end)) - Tstamp(ploc0(end)));
                    peaks_lag(dataIndex) = Tstamp(ploc1(end)) - Tstamp(ploc0(end));
%                     total_lag_peaks_arr(dataIndex) = tot_lag_peaks;
                    %dt_trough_tot_lag = dt_trough_tot_lag + (Tstamp1(ndploc1(end)) - Tstamp(ndploc0(end)));
                    locIndex = min(length(ploc1),length(ploc0));
                end


            end
%             avg_width0 = (avg_width0+w0(locIndex))/locIndex;
            prev_data0_peaks = curr_num_peaks0;
            prev_data0_peaks
            prev_data1_peaks = curr_num_peaks1;
%             meanCycle0 = mean(diff(ploc0));
%             meanCycle1 = mean(diff(ploc1));
         else
            if dataIndex>=2
%                 total_lag_peaks_arr(dataIndex) = total_lag_peaks_arr(dataIndex-1);
            end
        end
   
    end

%     Finds the time lag calculated from the peaks of the derivative
    if ~isempty(dploc0)&&~isempty(dploc1)
        if (curr_num_dpeaks1>prev_data1_dpeaks)||(curr_num_dpeaks0>prev_data0_dpeaks)
            if (length(dploc0)&&length(dploc1))==dlocIndex
%                 tot_lag_peaks = tot_lag_peaks + (ploc1(end) - ploc0(end) + (dploc1(end) - dploc0(end)))/2;
                dt_peak_tot_lag = dt_peak_tot_lag + (Tstamp(dploc1(dlocIndex)) - Tstamp(dploc0(dlocIndex)));
                dt_peaks_lag(dataIndex) = Tstamp(dploc1(dlocIndex)) - Tstamp(dploc0(dlocIndex));
                dlocIndex = dlocIndex + 1;
            else
                dt_peak_tot_lag = dt_peak_tot_lag + (Tstamp(dploc1(end)) - Tstamp(dploc0(end)));
                dt_peaks_lag(dataIndex) = Tstamp(dploc1(end)) - Tstamp(dploc0(end));
                %dt_trough_tot_lag = dt_trough_tot_lag + (Tstamp1(ndploc1(end)) - Tstamp(ndploc0(end)));
                dlocIndex = min(length(dploc1),length(dploc0));


            end
            prev_data0_dpeaks = curr_num_dpeaks0;
            prev_data1_dpeaks = curr_num_dpeaks1;
        end
    end
   
    addpoints(line0,Tstamp(dataIndex),data0);
    
    addpoints(line1,Tstamp(dataIndex),data1);
%     if abs(peaks_lag(dataIndex))>=1
%         color_tracer(dataIndex) = dataIndex;
%         addpoints(line2r,Tstamp(color_tracer(dataIndex)),peaks_lag(color_tracer(dataIndex)));
%         clearpoints(line2);
%         line2.Color = 'red';
% %         addpoints(line2r,Tstamp(dataIndex),peaks_lag(dataIndex));
%     elseif (1>abs(peaks_lag(dataIndex)))&&(abs(peaks_lag(dataIndex))>=0.5)
% %         clearpoints(line2)
%         line2.Color = 'magenta';
% %         addpoints(line2r,Tstamp(dataIndex),peaks_lag(dataIndex));
%     elseif (0.5>abs(peaks_lag(dataIndex)))&&(abs(peaks_lag(dataIndex))>=0.250)
% %         clearpoints(line2)
%         line2.Color = 'yellow';
% %         addpoints(line2r,Tstamp(dataIndex),peaks_lag(dataIndex));
%     else
% %         clearpoints(line2)
%         line2.Color = 'green';
% %         addpoints(line2r,Tstamp(dataIndex),peaks_lag(dataIndex));
%     end
    % Change color of time difference plot
    if abs(peaks_lag(dataIndex))>=1
        line2.Color = 'red';
    elseif (1>abs(peaks_lag(dataIndex)))&&(abs(peaks_lag(dataIndex))>=0.5)
        line2.Color = 'magenta';
    elseif (0.5>abs(peaks_lag(dataIndex)))&&(abs(peaks_lag(dataIndex))>=0.250)
        line2.Color = 'yellow';
    else
        line2.Color = 'green';
    end
    % Add peak points to live pressure sig0 and sig1
    if ~isempty(ploc0)
        clearpoints(line0p)
        addpoints(line0p,Tstamp(ploc0),peaks0);
    end
    if ~isempty(ploc1)
        clearpoints(line1p)
        addpoints(line1p,Tstamp(ploc1),peaks1);
    end

    
    addpoints(line2,Tstamp(dataIndex),peaks_lag(dataIndex));
    xlim([Tstamp(dataIndex)-dx, Tstamp(dataIndex)+(dx/2)]);

    drawnow %limitrate %nocallbacks

    dataIndex = dataIndex + 1;
    loop_counter= loop_counter+1;
end
dataIndex
Tstamp = nonzeros(Tstamp);
% total_lag_peaks_arr(end) = [];%total_lag_peaks_arr(end-1);
% peaks_lag(end) = peaks_lag(end);
% length_of_peaks_lag_arr = length(peaks_lag);
% total_lag_peaks_arr(end) = total_lag_peaks_arr(end-1);%total_lag_peaks_arr(end-1);
peaks_lag = nonzeros(peaks_lag);
% total_lag_peaks_arr = nonzeros(total_lag_peaks_arr);
avg_of_peaks_arr = mean(peaks_lag);
max_lag = max(peaks_lag); % Maximum instantaneous lag in seconds
min_lag = min(peaks_lag);
sum_peaks_lag = sum(peaks_lag); % total Time lag in seconds should == peaks_tot_lag

% freq0 = length(peaks0)/timer % strokes (peaks) per second
% freq1 = length(peaks1)/timer

% ppos0 = length(ploc0); % total number of peaks
% ppos1 = length(ploc1);


% Show peaks on plot

figure;
plot(Tstamp,deriv0(1:length(Tstamp)),"Color",'red');
hold on
plot(Tstamp,A0Array(1:length(Tstamp)),"Color",'blue');
plot(Tstamp,A1Array(1:length(Tstamp)),"Color",'cyan')
plot(Tstamp,deriv1(1:length(Tstamp)),"Color",'blue')
plot(Tstamp(ploc0),'go')
plot(Tstamp(ploc1),'mo')
plot(Tstamp(dploc0),dpeaks0,'g^')
plot(Tstamp(dploc1),dpeaks1,'m^')
plot(Tstamp(ndploc0),-ndpeaks0,'gv')
plot(Tstamp(ndploc1),-ndpeaks1,'rv')


% % Plot of curves with profile markings
figure; 
plot(Tstamp,A0Array(1:length(Tstamp)),"Color",'red'); hold on; plot(Tstamp(ploc0),peaks0,'go');
plot(Tstamp,A1Array(1:length(Tstamp)),"Color",'blue'); plot(Tstamp(ploc1),peaks1,'mo');
plot(Tstamp(dploc0),dpeaks0,'g^'); plot(Tstamp(dploc1),dpeaks1,'m^'); 
plot(Tstamp(ndploc0),ndpeaks0,'gv'); plot(Tstamp(ndploc1),ndpeaks1,'mv'); 
plot(Tstamp(start0),A0Array(start0),'g>',Tstamp(end0),A0Array(end0),'g<');
plot(Tstamp(start1),A1Array(start1),'m>',Tstamp(end1),A1Array(end1),'m<');

figure; 
plot(Tstamp,total_lag_peaks_arr(1:length(Tstamp)),"Color",'red')



% Find the sync factor, should add a length check
% for k = 1:length(ploc0)
%     sync(k) = (Tstamp1(ploc1(k)) - Tstamp(ploc0(k))) + (Tstamp1(dploc1(k)) - Tstamp(dploc0(k))) + (Tstamp1(ndploc1(k)) - Tstamp(ndploc0(k)));
% end

% print
avg_of_peaks_arr
max_lag
min_lag
sum_peaks_lag
tot_lag_peaks
dt_peak_tot_lag
% dt_trough_tot_lag
idx = find(peaks_lag==min_lag);
time_of_min = (idx);
time_of_min
% sum(sync)/3
smallest_array = min(length(ploc1),length(ploc0))
mallest_array2 = min(length(dploc1),length(dploc0))
smallest_array3 = min(length(ndploc1),length(ndploc0))
(length(ploc1)==length(ploc0))&&(length(dploc1)==length(dploc0))&&(length(ndploc1)==length(ndploc0))
for j = 1:smallest_array-1
    mid_0curve(j) = (Tstamp(dploc0(j))-Tstamp(ndploc0(j)))/2 + Tstamp(dploc0(j));
    mid_1curve(j) = (Tstamp(dploc1(j))-Tstamp(ndploc1(j)))/2 + Tstamp(dploc1(j));
end

% clear all;
%SADC(abs(SADC-mean(SADC))<2*std(SADC)) = 0; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
