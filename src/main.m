clear board0;
% Initialize the Arduino board

% board0 = serialport("COM5",19200);
% board1 = serialport('COM5',115200);

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
% % Loop control
duration = 10; % Time of demonstration
stepTime = 0.0001; % Time interval between consecutive data reads.
% samples = duration/stepTime; % Total number of data samples to be recorded.
% 
prev_0strokes = 0;
prev_1strokes = 0;
% prev_0dstrokes = 0;
% prev_1dstrokes = 0;
% 
% % Indices
dataIndex = 1;
locIndex = 1;
dlocIndex = 1;
% 
% % Array Definitions
% t = zeros(1,10000); % Array to track each looped peaks time difference
% nt = zeros(1,10000); % Array to track each looped trough time difference
% dt = zeros(1,10000); % Array to track each looped derivitive time difference
A0Array = zeros(1,200*duration);
A1Array = zeros(1,200*duration);
Tstamp0 = zeros(1,200*duration);
% % Tstamp1 = zeros(1,100000);
% 
% Initialize live plot
% start = datenum(datetime('now','Format','HH:mm:ss.SSS'));

% subplot(3, 1, 1)
line0 = animatedline('Color','r');
line1 = animatedline('Color','b');
% hold on
% title("Recorded Data");
% xlabel("Time (sec)");

% legend('Signal 0', 'Signal 1');
% ylim([0,5]);
% subplot(3, 1, 2)
% line2 = line(start-start,tlag,'Color','green');% Used to be time_diff
% % ylim([-2,2]);
% title("Time Lag");
% xlabel("Time");
% subplot(3, 1, 3)
% line3 = line(start-start,dtlag,'Color','magenta');
% title("Time Lag from Derivative");

tObj = tic; %Start Timer
while (toc(tObj)<= duration)

%     Tstamp1(dataIndex) = toc(tObj);%(datenum(datetime('now','Format','ss.SSS'))-start)*24*3600;%, 'Format', 'HH:mm:ss.SSS');
    data = fscanf(board0,'%d%d%d');
    time = data(1);
    data0 = data(2);
    data1 = data(3);
    A0Array(dataIndex) = data0;
    Tstamp0(dataIndex) = time;
    A1Array(dataIndex) = data1;


    deriv0 = gradient(A0Array);
    deriv1 = gradient(A1Array);
%     
%     % Find edges of curves
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

%      %Find Peaks
    [peaks0, ploc0] = findpeaks(A0Array,"MinPeakDistance",60);%,"MinPeakDistance",6,'MinPeakHeight',max(A0Array)-std(A0Array),'MinPeakProminence',3);%;%max(A0Array)-2*std(A0Array),'WidthReference','halfheight');
    [peaks1, ploc1] = findpeaks(A1Array,"MinPeakDistance",60);%,"MinPeakDistance",6,'MinPeakHeight',max(A1Array)-std(A0Array),'MinPeakProminence',3);%2*std(A0Array),'WidthReference','halfheight');
%     [dpeaks0, dploc0,dw0,dp0] = findpeaks(deriv0,"MinPeakDistance",6,'MinPeakHeight',max(deriv0)-3*std(deriv0));%,'MinPeakProminence',1.5);
%     [dpeaks1, dploc1,dw1,dp1] = findpeaks(deriv1,"MinPeakDistance",6,'MinPeakHeight',max(deriv1)-3*std(deriv1));%,'MinPeakProminence',1.5);
%     % Troughs
% %     [npeaks0, nploc0] = findpeaks(-A0Array,"MinPeakDistance",4,'MinPeakHeight',max(-A0Array)-2*std(-A0Array),'MinPeakProminence',0.1);
% %     [npeaks1, nploc1] = findpeaks(-A1Array,"MinPeakDistance",4,'MinPeakHeight',max(-A1Array)-2*std(-A1Array),'MinPeakProminence',0.1);
%     [ndpeaks0, ndploc0,ndw0,ndp0] = findpeaks(-deriv0,"MinPeakDistance",6,'MinPeakHeight',max(-deriv0)-3*std(-deriv0),'MinPeakProminence',1);%,'MinPeakProminence',0.1);
%     [ndpeaks1, ndploc1,ndw1,ndp1] = findpeaks(-deriv1,"MinPeakDistance",6,'MinPeakHeight',max(-deriv1)-3*std(-deriv1),'MinPeakProminence',1);%,'MinPeakProminence',0.1);

    num_strokes0 = length(ploc0);
    num_strokes1 = length(ploc1);
%     dnum_strokes0 = length(dploc0);
%     dnum_strokes1 = length(dploc1);
% 
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
% %             if ~(isempty(nploc1)&&isempty(nploc0))
% %                 ntlag = ntlag + (Tstamp1(nploc1(end)) - Tstamp0(nploc0(end)));
% %                 nt(dataIndex) = Tstamp1(nploc1(end))-Tstamp0(nploc0(end));
% %             end
% 
%             if ~isempty(dploc0)&&~isempty(dploc1)%length(ploc0)>=2
%                 if (dnum_strokes1>prev_1dstrokes)||(dnum_strokes0>prev_0dstrokes)
%                     if (length(dploc0)||length(dploc1))==dlocIndex
%     %                 tlag = tlag + (ploc1(end) - ploc0(end) + (dploc1(end) - dploc0(end)))/2;
%                         dtlag = dtlag + (Tstamp0(dploc1(dlocIndex)) - Tstamp0(dploc0(dlocIndex)));
%                         dt(dataIndex) = Tstamp0(dploc1(dlocIndex)) - Tstamp0(dploc0(dlocIndex));
%                         dlocIndex = dlocIndex + 1;
%                     else
%                         dtlag = dtlag + (Tstamp0(dploc1(end)) - Tstamp0(dploc0(end)));
%                         dt(dataIndex) = Tstamp0(dploc1(end)) - Tstamp0(dploc0(end));
%                         %ndtlag = ndtlag + (Tstamp1(ndploc1(end)) - Tstamp0(ndploc0(end)));
%                         dlocIndex = max(length(dploc1),length(dploc0));
%                     prev_0dstrokes = dnum_strokes0;
%                     prev_1dstrokes = dnum_strokes1;
%     
%                     end
%                 end
%             end 
        end
    end

    pause(stepTime);
    addpoints(line0,time,data0);
    addpoints(line1,time,data1);
%     subplot(3, 1, 1);


%     subplot(3, 1, 2);
%     line2.XData = [line2.XData time];
%     line2.YData = [line2.YData tlag];
% %     line2.YData = [line2.YData t(dataIndex)];
%     subplot(3, 1, 3);
%     line3.XData = [line3.XData time];
%     line3.YData = [line3.YData dtlag]; % Lag calc from Derivative
% %     line3.YData = [line3.YData data0-data1]; %Pressure difference, might be more useful


%     ARD_DATA_1 = timetable(A0Array',A1Array','TimeStep',seconds(stepTime)); This Works
%     plot(ARD_DATA_1.Time, ARD_DATA_1.Variables);
    drawnow;
    dataIndex = dataIndex + 1;
    

end
% t = nonzeros(t);
% avg = mean(t);
% maxi = max(t); % Maximum instantaneous lag in seconds
% mini = min(t);
% sumt = sum(ndtlag); % total Time lag in seconds
% 
% freq0 = length(peaks0)/timer % Number of strokes
% freq1 = length(peaks1)/timer
% ppos0 = length(ploc0);
% ppos1 = length(ploc1);

% Show peaks on plot
% subplot(3,1,3)
% EndPlots = figure(2);
% ax = axes('Parent', EndPlots);
% % plot(A0Array,"Color",'blue');
% plot(Tstamp0,deriv0,"Color",'red');
% hold on
% % plot(A1Array,"Color",'cyan')
% plot(Tstamp1,deriv1,"Color",'blue')
% % plot(ploc0,peaks0,'r^')
% % plot(ploc1,peaks1,'g^')
% plot(dploc0,dpeaks0,'g^')
% plot(dploc1,dpeaks1,'r^')
% plot(ndploc0,-ndpeaks0,'g^')
% plot(ndploc1,-ndpeaks1,'r^')
% 
% % Plot of curves with profile markings
figure; 
plot(Tstamp0,A0Array,"Color",'red'); hold on; plot(Tstamp0(ploc0),peaks0,'go');
plot(Tstamp0,A1Array,"Color",'blue'); plot(Tstamp0(ploc1),peaks1,'mo');
% plot(Tstamp0(dploc0),dpeaks0,'g^'); plot(Tstamp0(dploc1),dpeaks1,'m^'); 
% plot(Tstamp0(ndploc0),ndpeaks0,'gv'); plot(Tstamp1(ndploc1),ndpeaks1,'mv'); 
plot(Tstamp0(pstart0),A0Array(pstart0),'g>',Tstamp0(pend0),A0Array(pend0),'g<');
plot(Tstamp0(pstart1),A1Array(pstart1),'m>',Tstamp0(pend1),A1Array(pend1),'m<');
% 
% for k = 1:length(ploc0)
%     sync(k) = (Tstamp1(ploc1(k)) - Tstamp0(ploc0(k))) + (Tstamp1(dploc1(k)) - Tstamp0(dploc0(k))) + (Tstamp1(ndploc1(k)) - Tstamp0(ndploc0(k)));
% end
% avg
% maxi
% mini
% sumt
% tlag
% dtlag
% ndtlag
% idx = find(t==mini);
% time_of_min = (idx);
% time_of_min
% sum(sync)/3
% 
% hold off



%SADC(abs(SADC-mean(SADC))<2*std(SADC)) = 0; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
