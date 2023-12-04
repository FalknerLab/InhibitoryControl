%created 10-18-19 %modified 03-10-20 %author e mae guthman

%edit analysis of mean responses to mean sweep rather than aligned by max rise

close all
clearvars

%% INIT VARS
%mpoa or shell input (1 mpoa, 0 sh)
choice = 1;

%sample rate
fs = 10000; %Hz

%stimulus info in ms
stimstart =  925*(fs/1000); %if light duration 2ms; 5638.6*(fs/1000); %if light duration 15ms 
stimdur = 1*(fs/1000); %if light duration 2ms; 15*(fs/1000); %if light duration 15ms 
stimend = stimstart+stimdur;

%set baseline -501ms to -1ms before stim
blstart=stimstart-500*(fs/1000)-1*(fs/1000);
blend=stimstart-(fs/1000);

%plot color
if choice == 1
    cc = [180 151 214]./255; %wisteria
elseif choice == 0
    cc = .8.*[89 217 153]./255 %medium aquamarine
end
chrColor = [85 155 250]./255; %blue

%% LOAD DATA
[file, path] = uigetfile('*.abf');
cd(path);
oStimFile = fullfile(path,file);

rawData = abfload(oStimFile,'sweeps','a');
[b, a] = butter(1, 2000/fs, 'low'); %2 kHz lowpass
filtData(:,:) = filtfilt(b, a, rawData(:,1,:));

%time vector
t = 1/(fs/1000):1/(fs/1000):size(filtData,1)/(fs/1000);

%% ANALYSIS
% remove Na-spike sweeps/Ra >35
% 65sh_3: [2 3 4 5 6 10 17 18 20];
% 66sh_1: [2 4 5 10 11 12 13 15 17 18 19 20 22 23 25 28 29 30]
% 93mpoa2: [5 6 8 9 11 14 15 16 17 20]
% 93mpoa4: [2 4 7 9 11 13 16]
% 17sh4: [9]
% 13mpoa3: [2 12 18]
% 78mpoa1: [6 8 14 31] 
% 78mpoa3: [5 14 16 18 19]
% 62mpoa1: [1]
% 61mpoa1: [2 6 10 11 12 13]
% 
rawData(:,:,[1]) = []; 
filtData(:,[1]) = [];

%CHECK RA
%injection for Ra
injstart = 225.2*(fs/1000); %voltage injection start
injend = 275.1*(fs/1000); %voltage injection end

%Ra analysis
[Ra, meanRa, keep] = getRaThirty(rawData,injstart,injend,fs)  %#ok<*NOPTS>

% if keep == 1
    %baseline subtract traces
    blData = -1*ones(size(filtData));
    for ii = 1:size(filtData,2)
        blData(:,ii) = filtData(:,ii) - mean(filtData(blstart:blend,ii));
    end
    
    %get properties of events
    [output.amplitude, output.latency, output.rise, output.decay, output.tAllPeaks] = getOptStimCurrent(blData,stimstart,stimdur,fs);
    
    if sum(isnan(output.amplitude)) ~= length(output.amplitude)
        %mean amplitude
        amp=output.amplitude;
        amp(isnan(amp))=[];
        output.mAmp=mean(amp);
        
        %latency info
        lat=output.latency;
        lat(isnan(lat))=[];
        output.mLat=mean(lat);
        output.jitter=std(lat);
        
        %kinetics
        rt = output.rise;
        rt(isnan(rt)) = [];
        output.mRT = mean(rt);
        output.decay;
        
        %success rate
        noFail=sum(isnan(output.amplitude));
        output.sRate=100*((length(output.amplitude)-noFail)/length(output.amplitude));
        
        %mean trace
        output.mTrace = mean(blData(:,~isnan(output.amplitude)),2);
        output.alignedTraces = blData(:,~isnan(output.amplitude));
        
        %get values of failures
        windowCenter = mean(output.tAllPeaks(~isnan(output.tAllPeaks)));
        ampFailures = zeros(1,length(output.amplitude));
        for ii = 1:length(output.amplitude)
            if isnan(output.amplitude(ii))
                ampFailures(ii)=min(blData(windowCenter-.005*fs:windowCenter+.005*fs,ii));
            end
        end
    else
        %success rate
        output.sRate=0;
    end
% end

%% PLOT
oStimFig = figure(1);
oStimFig.Position = [20 810 1400 475];

%current traces
subplot(1,9,1:7)
hold on
if sum(isnan(output.amplitude))>0
    plot(t,blData(:,isnan(output.amplitude(:))),'color','k','linewidth',1)
end
for ii = 1:size(blData,2)-noFail
    plot(t(1:length(output.mTrace)),output.alignedTraces(:,ii),'color',cc./1.5,'linewidth',.5)
end
plot(t(1:length(output.mTrace)),output.mTrace,'color',cc,'linewidth',3)
line([stimstart/(fs/1000) stimstart/(fs/1000)+stimdur/(fs/1000) ], [output.mAmp/-10 output.mAmp/-10],'color',chrColor,'linewidth',10)
xlim([stimstart/(fs/1000)-5 stimstart/(fs/1000)+70])
ylabel('current (pA)')
xlabel('time (ms)')
oStimAx = gca;
setAx(oStimAx);

%amplitudes
subplot(1,9,8:9)
hold on
for ii = 1:size(filtData,2)
    if isnan(output.amplitude(ii))
        scatter(ii,ampFailures(ii),150,'k','filled')
    else
        scatter(ii,output.amplitude(ii),150,cc,'filled')
    end
end
ylim(oStimAx.YLim)
xlabel('trial')
setAx(gca);

output.vector =[output.mAmp, output.mLat, output.jitter output.mRT, output.decay]';
output.labels = {'oIPSC amplitude', 'oIPSC latency', 'oIPSC jitter', 'oIPSC 20-80% RT', 'oIPSC decay'};

%% SAVE
if choice == 0
    svDir = '/Volumes/falkner/Mae/patch_data/2021_mpoa_vmhshell_stim/analyzed_data/oStim/single_pulse/shell_stim/'
elseif choice == 1
    svDir = '/Volumes/falkner/Mae/patch_data/2021_mpoa_vmhshell_stim/analyzed_data/oStim/single_pulse/mpoa_stim/'
end
appendThis = input('file name for save (in quotes):');
save([svDir, appendThis,  'analyzed.mat'],'output','oStimFig')