%created 08-12-2021 %author eartha mae guthman
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
    cc = .9.*[89 217 153]./255 %medium aquamarine
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
% remove Na-spikes + Ra outside +/-20%
% 93mpoa5: [6 9 11 15 20]
% 21sh3: [4]
% 17sh1: [4]
% 17sh4: [3 5 8 11]
% 13mpoa3: [20]
% 78mpoa1: [4 12 18 19 21 25 27]
% 78mpoa3: [11]
% 67sh2: [8 9]
% 62mpoa1: [4 9 10];
% 61mpoa1: [10]
% 24sh2: [10 11 13 15 17]
% 26sh1: [5 7 9 10]
% 26sh2: [1 2 19]

% rawData(:,:,[4 9 10]) = []; 
% filtData(:,[4 9 10]) = [];

%CHECK RA
%injection for Ra
injstart = 225.2*(fs/1000); %voltage injection start
injend = 275.1*(fs/1000); %voltage injection end

%Ra analysis
[Ra, meanRa, keep] = getRaThirty(rawData,injstart,injend,fs)  %#ok<*NOPTS>

if keep == 1
    %baseline subtract traces
    blData = -1*ones(size(filtData));
    for ii = 1:size(filtData,2)
        blData(:,ii) = filtData(:,ii) - mean(filtData(blstart:blend,ii));
    end
    
    %get mean trace
    meanData = mean(blData,2);
    
    % get PSC peaks
    %set windows for capturing peak
    dpISI = fs*0.02; %isi in data ponts
    msISI = 20; %isi ms
    windowS = 0.001*fs; %in dp start window 1ms after stimulus
    windowE = dpISI-0.001*fs; %window for peak is 2ms post stimulus to 1ms pre next stimulus
    
    %find peaks
    for jj = 1:2
        [peakVals(jj), peakInds(jj)] = min(meanData(stimstart+windowS+dpISI*(jj-1):stimstart+windowE+dpISI*(jj-1)));
        peakInds(jj) = peakInds(jj) + (stimstart+windowS+dpISI*(jj-1)) - 1; %fix indexing
    end
    
    % get paired pulse ratio
    output.ppr = peakVals(2)/peakVals(1);
end

%% PLOT
ppStimFig = figure(1);
ppStimFig.Position = [20 810 375 475];

%current traces
hold on
for ii = 1:size(blData,2)
    plot(t,blData(:,ii),'color',cc./1.5,'linewidth',.5)
end
plot(t,meanData,'color',cc,'linewidth',3)
line([stimstart/10 stimstart/10+stimdur/10], [50 50],'color',chrColor,'linewidth',10)
line([stimstart/10+20 stimstart/10+20+stimdur/10], [50 50],'color',chrColor,'linewidth',10)
xlim([stimstart/10-3 stimstart/10+50])
ylabel('current (pA)')
xlabel('time (ms)')
ppStimAx = gca;
setAx(ppStimAx);

%% SAVE
if choice == 0
    svDir = '/Volumes/falkner/Mae/patch_data/2021_mpoa_vmhshell_stim/analyzed_data/oStim/paired_pulse/shell_stim/'
elseif choice == 1
    svDir = '/Volumes/falkner/Mae/patch_data/2021_mpoa_vmhshell_stim/analyzed_data/oStim/paired_pulse/mpoa_stim/'
end
appendThis = input('file name for save (in quotes):');
save([svDir, appendThis,  '_paired_pulse_analyzed.mat'])