clc;
close all;
clear board0 dataFromCSV;
% Initialize the Arduino board

useArduino = true;
if useArduino
    board0 = serialport("COM3", 19200);
    configureTerminator(board0, "CR/LF");
    % board1 = serialport('COM5',115200);
else
    tableFromCSV = readtable(['..\data\data_5_modified.csv']);
    dataFromCSV = table2array(tableFromCSV);
end

data0 = 0;
data1 = 0;
time = 0;

% Scalars
tlag = 0; % Total lag gotten from peaks
ntlag = 0; % negative waveform, Total lag from troughs
dtlag = 0; % Derivative, Total lag from derivitive peaks
ndtlag = 0; % troughs of derivative
sync = 0;
tol = 0.05;

%% Loop control
duration = 15; % Time of demonstration
stepTime = 0.0001; % Time interval between consecutive data reads.

prev_0strokes = 0;
prev_1strokes = 0;

%% Indices
dataIndex = 1;
locIndex = 1;
dlocIndex = 1;

%% Array Definitions
A0Array = zeros(1,200*duration);
A1Array = zeros(1,200*duration);
Tstamp0 = zeros(1,200*duration);

line0 = animatedline('Color','r');
line1 = animatedline('Color','b');

tObj = tic; % Start Timer
index = 0; % used for CSV input only
ticStart = tic;
while (toc(tObj)<= duration)
tic
if useArduino
    raw = writeread(board0, "a");
    data = sscanf(raw, '%d');
else
    index = index + 1;
    data = dataFromCSV(index, :)
end

    time = data(1);
    data0 = data(2);
    data1 = data(3);
    A0Array(dataIndex) = data0;
    Tstamp0(dataIndex) = time;
    A1Array(dataIndex) = data1;


    deriv0 = gradient(A0Array);
    deriv1 = gradient(A1Array);

%% Find edges of curves
    pbool0L = (tol*max(A0Array)<A0Array);
    pbool0G = (A0Array>min(A0Array));
    pbool0 = pbool0L.*pbool0G;
    pbool0 = [false,pbool0,false];
    pstart0 = strfind(pbool0,[0 1]);
    pend0 = strfind(pbool0,[1 0]);
    pbool1L = (tol*max(A1Array)<A1Array);
    pbool1G = (A1Array>min(A1Array));
    pbool1 = pbool1L.*pbool1G;
    pbool1 = [false,pbool1,false];
    pstart1 = strfind(pbool1,[0 1]);
    pstart0 = pstart0 - 1;
    pstart1 = pstart1 - 1;
    pend1 = strfind(pbool1,[1 0]);

%% Find Peaks
    [peaks0, ploc0] = findpeaks(A0Array,"MinPeakDistance",60);%,"MinPeakDistance",6,'MinPeakHeight',max(A0Array)-std(A0Array),'MinPeakProminence',3);%;%max(A0Array)-2*std(A0Array),'WidthReference','halfheight');
    [peaks1, ploc1] = findpeaks(A1Array,"MinPeakDistance",60);%,"MinPeakDistance",6,'MinPeakHeight',max(A1Array)-std(A0Array),'MinPeakProminence',3);%2*std(A0Array),'WidthReference','halfheight');

    num_strokes0 = length(ploc0);
    num_strokes1 = length(ploc1);

    if ~isempty(ploc0)&&~isempty(ploc1)
        if (num_strokes1>prev_1strokes)||(num_strokes0>prev_0strokes)
            if length(ploc0)==length(ploc1)
                tlag = tlag + (Tstamp0(ploc1(locIndex)) - Tstamp0(ploc0(locIndex)));
                t(dataIndex) = Tstamp0(ploc1(locIndex))-Tstamp0(ploc0(locIndex));
                locIndex = locIndex + 1;
            else
                tlag = tlag + (Tstamp0(ploc1(end)) - Tstamp0(ploc0(end)));
                locIndex = max(length(ploc1),length(ploc0));
            end

            prev_0strokes = num_strokes0;
            prev_1strokes = num_strokes1;
        end
    end

    pause(stepTime);
    addpoints(line0,time,data0);
    addpoints(line1,time,data1);

    drawnow;
    dataIndex = dataIndex + 1;
    pause(0.04);
    %toc
end

%% Plot of curves with profile markings
figure; 
plot(Tstamp0,A0Array,"Color",'red'); 
hold on; 
plot(Tstamp0(ploc0),peaks0,'go');
plot(Tstamp0,A1Array,"Color",'blue'); 
plot(Tstamp0(ploc1),peaks1,'mo');

plot(Tstamp0(pstart0),A0Array(pstart0),'g>',Tstamp0(pend0),A0Array(pend0),'g<');
plot(Tstamp0(pstart1),A1Array(pstart1),'m>',Tstamp0(pend1),A1Array(pend1),'m<');
