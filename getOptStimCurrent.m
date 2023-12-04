function [amplitude, latency, rise, decay, tPeaks] = getOptStimCurrent(data,stimstart,stimdur,samplerate)
% Function to get data on optically evoked post-synaptic currents, created 10-18-2019,
% modified 10-03-2021
% takes sweeps and outputs data on the evoked post-synaptic currents
% returns PSC amplitude, latency, 20-80 rise time, weighted decay tau

%% inits
%input inits
if length(size(data)) == 2
    noSweeps=size(data,2);
elseif length(size(data)) > 2
    error('Too many input arguments in data file')
end

%baseline inits -- 500ms before 1ms before stimulus start
blstart=stimstart-0.5*samplerate-1*(samplerate/1000);
blend=blstart+0.5*samplerate;

%stim inits
stimend=stimstart+stimdur;

%create vectors for output
amplitude=zeros(noSweeps,1);
latency=zeros(noSweeps,1);
rawEventStart=zeros(noSweeps,1);
rise=zeros(noSweeps,1);

%create other vectors
peakRawCurrent=zeros(noSweeps,1);
tPeaks=zeros(noSweeps,1);

%% find data
%get data on the currents
for ii = 1:noSweeps
    %amplitude, 10ms window
    [amplitude(ii),peakRawCurrent(ii),tPeaks(ii)] = findPSCpeak(data,stimstart,.01*samplerate+stimdur,blstart,blend,ii);
    amplitude(amplitude==0) = NaN;
    
    %get rise time and latency
    %rise time: 20-80% rise time
    %latency: time from stimulus end to 20% rise
    if isnan(amplitude(ii)) ~= 1 %if event occurred
        dataForRiseLat = data(stimend:tPeaks(ii),ii); %get window of data from light off to peak
        dataForRiseLat = dataForRiseLat - mean(data(blstart:blend,ii)); %baseline subtract data
        [rise(ii), ~, latency(ii)] = findRise(dataForRiseLat,1,length(dataForRiseLat),samplerate);
        latency(ii) = latency(ii) - 1; %scales back to end of stimulus
        latency(ii) = latency(ii)*(1000/samplerate); %converts to ms
    elseif isnan(amplitude(ii)) %if no event occurred
        latency(ii) = NaN;
        rise(ii) = NaN;
    end
    
    
    %check that correct peak is found
    viewsweepsFig=figure('Position',[20 700 900 600]);
    hold on
    plot(data(:,ii))
    scatter(tPeaks(ii),peakRawCurrent(ii))
    if isnan(latency(ii)) ~= 1
        scatter(latency(ii)*(samplerate/1000)+stimstart-1,data(latency(ii)*(samplerate/1000)+stimstart-1),ii)
    end
    line([stimstart-.002*samplerate stimend+.05*samplerate],...
        [mean(data(blstart:blend,ii))-5*getMAD(data(blstart:blend,ii)) ...
        mean(data(blstart:blend,ii))-5*getMAD(data(blstart:blend,ii))],'color','k')
    xlim([stimstart-.002*samplerate stimend+.05*samplerate])
    if isnan(amplitude(ii)) ~= 1 %if event occurred
        ylim([1.25*amplitude(ii) mean(data(blstart:blend,ii))+5*getMAD(data(blstart:blend,ii))])
    else
        ylim([-10*getMAD(data(blstart:blend,ii)) 10*getMAD(data(blstart:blend,ii))])
    end
    pfft = 1; %pause here to examine every sweep by eye for currents
    
    close(viewsweepsFig)
end

%get mean trace aligned max rise
if sum(~isnan(amplitude)) > 0 %if any successes get mean data and tau
    meanTrace = mean(data(:,~isnan(amplitude)),2);
    
    %weighted decay tau
    decay=findTau(meanTrace(1:stimstart+150*samplerate/1000),mean(latency(~isnan(latency)))*samplerate/1000+stimstart,samplerate,'Double');
else
    decay = NaN;
end

end

function [amp, peakRawCurrent, time, direction] = findPSCpeak(thisData,tStart,tWindow,blstart,blend,sweepNo)

%get shortened vector of data for analysis
blCurrent=mean(thisData(blstart:blend,sweepNo));
madBaseline=getMAD(thisData(blstart:blend,sweepNo)); %median absolute deviation instead of standard deviation as it isn't affected by large deflections from mean as much as sd is
data2analyze=thisData(tStart:tStart+tWindow,sweepNo)-blCurrent;

%find peaks
[valleys(:,1),valleys(:,2)]=findvalleys(data2analyze);
abovethresh=valleys(:,1)<-5*madBaseline; % find those valleys from above that are above thresh
if sum(abovethresh) > 0
    peaksAboveThresh(:,1)=valleys(abovethresh,1);
    peaksAboveThresh(:,2)=valleys(abovethresh,2);
    [amp,peakInd]=min(peaksAboveThresh(:,1));
    peakRawCurrent=amp+blCurrent;
    time=peaksAboveThresh(peakInd,2)+tStart-1;
    
    %check that correct peak is found
    viewsweepsFig=figure;
    hold on
    plot(data2analyze)
    scatter(peaksAboveThresh(peakInd,2),amp)
    line([1 length(data2analyze)],[-10*madBaseline -10*madBaseline],'color','k')
%     line([1 length(data2analyze)],[10*madBaseline 6*madBaseline],'color','k')
    close(viewsweepsFig)
else
    amp = NaN;
    peakRawCurrent=NaN;
    time = NaN;
end
end
