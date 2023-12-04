%created 10-18-19 %modified 08-02-21 %author eartha mae guthman
close all
clearvars

%% INIT VARS

saveOn = 1;

%sample rate
fs = 10000; %Hz

%stimulus start in ms
injstart = 547.0*(fs/1000);
injend = 1546.9*(fs/1000);

%baseline
blstart = injstart-500*(fs/1000)-1*(fs/1000);
blend = injstart-1*(fs/1000);

%plot color
plotColor = [101, 31, 112]./(255*2); %purple

%% LOAD DATA
[file, path] = uigetfile('*.abf');
cd(path)
ciFile = fullfile(path,file);


%% PREPARE DATA
output.ccFile = abfload(ciFile,'sweeps','a');
%creates 3D matrix (MxNxP) where M is time, N is input/command, P is sweep number
%M is dependent on sample rate
%first data point in pCFile is at t=0
noSweeps = 21; %size(output.ccFile,3); %gives number of sweeps per file
t = 1/fs:1/fs:size(output.ccFile,1)/fs;

%find dI
sweepDelta = zeros(noSweeps,2);
sweepDelta(:,1) = 1:1:noSweeps;
sweepDelta(:,2) = -200:20:200;
[~,sag_sweep] = min(abs(sweepDelta(:,2)+100)); %finds sweep index closest to -100 sweep

%find AP Threshold
[thVals,thTimes]=findSpkThresh(output.ccFile,injstart,injend,fs,sweepDelta(:,2),plotColor);
thTimes(thTimes~=-1) = thTimes(thTimes~=-1) +1;
rheoSweep=find(thVals(:,1)~=-1);
if length(rheoSweep) >1
    rheoSweep(2:end)=[];
    output.rheo = sweepDelta(rheoSweep,2);
elseif isempty(rheoSweep)
    rheoSweep=noSweeps+1;
end


%% Display Traces
ciFig = figure(1);
title('Voltage Responses to Current Injections')
ciFig.Position = [20 775 825 525];

%get color range
colorShift = 1-max(plotColor);
for ii = 1:3
    cRange(:,ii) = plotColor(ii):colorShift/noSweeps:plotColor(ii)+colorShift;
end

%plot
subplot(3,1,[1 2])
hold on
for ii = 1:noSweeps
    plot(t,output.ccFile(:,1,ii),'color',cRange(end-ii+1,:),'linewidth',1.5)
end
xlim([injstart./fs-.05 injend./fs+.1])
ylabel('Vm (mV)')
setAx(gca);
subplot(3,1,3)
hold on
for ii = 1:noSweeps
    plot(t,output.ccFile(:,2,ii),'color',cRange(end-ii+1,:),'linewidth',1.5)
end
xlim([injstart./fs-.05 injend./fs+.1])
xlabel('time (s)')
ylabel('injected current (pA)')
setAx(gca);

%% GET MEMBRANE PROPERTIES
%membrane resistance, get with low delta experiment if two runs
output.rm=getRm(output.ccFile,blstart,blend,injstart,fs,rheoSweep,sweepDelta(:,2),cRange);

%voltage sag and rebound spikes, taken at -180 pA sweep
output.vsag = getSag(output.ccFile,sag_sweep,injstart,injend,blstart,blend,fs,cRange);
[~, ~, output.nRbnd] = getRbndspks(output.ccFile,sag_sweep,injend+1,injend+.5*fs,cRange,fs);

%Max firing rate, defined as inverse of ISI during first 200ms of spiking in last current injection before reduction in AP firing occurs
[output.maxFR, maxSweep, output.isi_cv, isi_sweep, output.isi_one, output.isi_two, output.FR, output.delay] = getMaxFiring(output.ccFile,noSweeps,rheoSweep,thTimes,fs,injstart,cRange);

%F-I plot
output.FIcurve = getFICurve(output.ccFile,noSweeps,rheoSweep,sweepDelta,thTimes,fs,injstart,cRange);

%rheobase data
%set up figure
rheoFig = figure(7);
rheoFig.Position = [875 775 825 525];

%spike train
subplot(1,8,1:4)
hold on
plot(t,output.ccFile(:,1,isi_sweep),'linewidth',2,'color',cRange(isi_sweep,:))
setAx(gca)
xlim([injstart./fs-.05 injend./fs+.1])
xlabel('time (s)')
ylabel('Vm (mV)')
title('rheobase action potentials')

%spike overlays
subplot(1,8,5:6)
hold on
tRheoSpikes = thTimes(isi_sweep,thTimes(isi_sweep,:)>-1);
output.rWaves=findHypSpkWaves(output.ccFile,tRheoSpikes,isi_sweep,fs); %aligned to threshold
for ii = 1:size(output.rWaves,1)
    plot(1000/fs.*(1:size(output.rWaves,2)),output.rWaves(ii,:),'color',cRange(end-isi_sweep+1,:),'linewidth',2)
end
setAx(gca)
xlabel('time (ms)')
title('overlaid spikes')

%phase plot
subplot(1,8,7:8)
hold on
for ii = 1:size(output.rWaves,1)
output.dVdt(ii,:)=diff(output.rWaves(ii,:)).*10; %gives derivative in mV/ms
plot(output.rWaves(ii,2:end),output.dVdt(ii,:),'color',cRange(end-isi_sweep+1,:),'linewidth',2)
end
xlabel('Vm (mV)')
ylabel('dV/dt (mV/ms)')
title('phase plot')
phaseAx = gca;
setAx(phaseAx);
phaseAx.YAxisLocation= 'origin'; phaseAx.XAxisLocation = 'origin';

%find ap peak data for rheobase sweep
[output.valRheoPeaks, output.tRheoPeaks, output.noRheoSpks] = getPeaks(output.ccFile,thVals,thTimes,isi_sweep);

%latency to first AP; defined as time from current injection to first AP peak
output.spkLat = (output.tRheoPeaks(1) - injstart)*(1000/fs); %gives first spike latency in ms

%threshold at rheobase
output.rheoThreshold = thVals(isi_sweep,1:output.noRheoSpks);

subplot(1,8,5:6)
for ii = 1:output.noRheoSpks
    thesePts = (find(output.rWaves(ii,:)==output.rheoThreshold(ii)));
    scatter((thesePts(1)-1)*(1000/fs),output.rheoThreshold(ii),50,cRange(end,:),'filled')
    clear thesePts
end

subplot(1,8,7:8)
for ii = 1:output.noRheoSpks
    thesePts = (find(output.rWaves(ii,:)==output.rheoThreshold(ii)));
    scatter(output.rheoThreshold(ii),output.dVdt(ii,thesePts(1)-1),50,cRange(end,:),'filled')
    clear thesePts
end

%spk amplitude at rheobase
output.rheoAmp = output.valRheoPeaks-thVals(isi_sweep,1:output.noRheoSpks)';

%AHP
[output.AHP, output.tAHP, output.decayAHP] = getAHP(output.ccFile,output.tRheoPeaks,thTimes(isi_sweep,1:output.noRheoSpks),injend,isi_sweep,fs);

figure

output.rheoAHPlat=zeros(length(output.AHP),1);
output.rheoAHPval=zeros(length(output.AHP),1);
for ii = 1:length(output.AHP)
    output.rheoAHPlat(ii)=(output.tAHP(ii)-thTimes(isi_sweep,ii))/(fs/1000);
    output.rheoAHPval(ii)=thVals(isi_sweep,ii)-output.AHP(ii);
end

%halfwidth
output.halfwidth = getHalfwidth(output.ccFile,isi_sweep,fs,thTimes(isi_sweep,1:output.noRheoSpks),output.tRheoPeaks,output.rheoAmp);

%max FR data
%find ap peak data for maxFR sweep
[output.maxSweepPeakVals, output.maxSweepPeakTimes, output.maxSweepSpks] = getPeaks(output.ccFile,thVals,thTimes,maxSweep);

%FR adaptation, ratio of first first ISI to last 3 ISI (<1 slowing down; >1 speeding up)
[output.frRatio, ~] = getRatioFR(thTimes,maxSweep,output.maxSweepSpks);

%amplitude adaptation, ratio of last three AP amplitudes to first AP amplitude (<1 getting smaller; >1 getting bigger)
output.maxSweepAmps = output.maxSweepPeakVals-thVals(maxSweep,1:output.maxSweepSpks);
output.ampRatio= mean([output.maxSweepAmps(end) output.maxSweepAmps(end-1) output.maxSweepAmps(end-2)])/output.maxSweepAmps(1);

%AP broadening, defined in Keshavarzi et al 2014 as the change in hw from AP1 to AP2
%redefined here as the ratio between the two (<1 gets faster; >1 gets slower)
output.halfwidthMaxSw = getHalfwidth(output.ccFile,maxSweep,fs,thTimes(maxSweep,1:output.maxSweepSpks),output.maxSweepPeakTimes,output.maxSweepAmps);
output.broadeningRatio = output.halfwidthMaxSw(2)/output.halfwidthMaxSw(1);

%deltaAHP
[output.MaxSwAHP, output.tMaxSwAHP] = getAHP(output.ccFile,output.maxSweepPeakTimes,thTimes(maxSweep,1:output.maxSweepSpks),injend,maxSweep,fs);
output.dAHP = output.MaxSwAHP(end)-output.MaxSwAHP(1);

%post injection voltage step (mean membrane voltage change from the 150ms before the injection to the 150ms after)
output.v_step = mean(output.ccFile(injend:injend+0.15*fs,1,maxSweep))- mean(output.ccFile(injstart-0.15*fs:injstart,1,maxSweep));


%% colalte data output vector
output.labels = {'membrane_resistance', 'max_firing_rate','max_ap_slode','max_ap_downslope','first_ap_lat', 'ap_threshold', 'ahp_amplitude', 'ahp_latency','ahp_rise',...
    'ap_amplitude', 'ap_halfwidth', 'isi_coef_var', 'first_isi','second_isi','firing_bias','delta_ahp','post_injection_voltage_step','FR_accomodation_ratio','amp_accomodation_ratio','halfwidth_broadening_ratio'};
output.vector = [output.rm; output.maxFR; mean(max(output.dVdt')); mean(min(output.dVdt')); output.spkLat; mean(output.rheoThreshold); mean(output.rheoAHPval); mean(output.rheoAHPlat); output.decayAHP; ...
    mean(output.rheoAmp); mean(output.halfwidth); output.isi_cv; output.isi_one; output.isi_two; output.delay; output.dAHP; output.v_step; output.frRatio;  output.ampRatio; output.broadeningRatio];

output.table = table(output.vector,'RowNames',output.labels);

output.vector %#ok<*NOPTS>

%% save
if saveOn
    svDir = '/Volumes/falkner/Mae/patch_data/2021_mpoa_vmhshell_stim/analyzed_data/membrane_props/';
    appendThis = input('file name for save (in quotes):');
    save([svDir, appendThis,'_current_pulses'])
end

%% functions

function [thVals,thTimes] = findSpkThresh(file,injectionStart,injectionEnd,fs,dI,pColor)

%init vars
noSweeps=length(dI);
thVals=[];
thTimes=[];
dVm = diff(file); %takes approx derivative (lose element n = 1)
dVm(:,2,:) = [];
dVm = dVm.*(fs/1000); %converts to V/s
file(1,:,:) = []; %remove element n = 1 from Vm data
t = fs^-1:fs^-1:size(file,1)/fs;
dVmThresh = 20; %manually set to 20V/s

for ii = 1:noSweeps
    tThresh=[];
    oThresh=[];
    if dI(ii) > 0
        %find spikes
        noSpkSweep=1;
        oThresh=find(dVm((.002*fs+injectionStart):injectionEnd,1,ii)>dVmThresh); %find points above the threshold of dVm during current injection (2ms buffer from start of step)
        for jj = 1:length(oThresh)
            if jj == 1
                %get point before and after threshold
                p0 = oThresh(jj)+injectionStart+.002*fs-1; %after
                p1 = oThresh(jj)+injectionStart+.002*fs-2; %before
                afterVal = abs(dVmThresh-dVm(p0,1,ii));
                beforeVal = abs(dVmThresh-dVm(p1,1,ii));
                if afterVal <= beforeVal %if value after threshold is closer to threshold...
                    tThresh(noSpkSweep) = p0;
                else %if value before threshold is closer to threshold...
                    tThresh(noSpkSweep) = p1;
                end
                noSpkSweep=noSpkSweep+1;
            else
                if oThresh(jj)>.001*fs+tThresh(noSpkSweep-1)-injectionStart %if time over threshold is at least 1ms past the last AP (refractory period)
                    p0 = oThresh(jj)+injectionStart+.002*fs-1; %after
                    p1 = oThresh(jj)+injectionStart+.002*fs-2; %before
                    afterVal = abs(dVmThresh-dVm(p0,1,ii));
                    beforeVal = abs(dVmThresh-dVm(p1,1,ii));
                    if afterVal <= beforeVal %if value after threshold is closer to threshold...
                        tThresh(noSpkSweep) = p0;
                    else %if value before threshold is closer to threshold...
                        tThresh(noSpkSweep) = p1;
                    end
                    noSpkSweep=noSpkSweep+1;
                end
            end
        end
        noSpkSweep = noSpkSweep -1
        
        %Plots AP traces and thresholding
        figure(1)
        subplot(1,2,1)
        hold on
        plot(t,file(:,1,ii),'color',pColor.*1.5)
        ylim([min(min(file(:,1,:)))-0.1*max(range(file(:,1,:))) max(max(file(:,1,:)))+0.1*max(range(file(:,1,:)))])
        ylabel('Membrane Potential (mV)')
        xlabel('Time (ms)')
        xlim([.5 1.75])
        if isempty(tThresh) ~= 1
            scatter(tThresh/fs,file(tThresh,1,ii),'markeredgecolor',pColor)
        end
        legend('Vm','threshold');
        thisLegend = legend;
        thisLegend.Box = 'off';
        thisAx = gca;
        setAx(thisAx);
        subplot(1,2,2)
        hold on
        plot(t,dVm(:,1,ii),'color',pColor.*1.5)
        line([t(1) t(end)],[dVmThresh dVmThresh],'color','k','linewidth',1.25);
        ylabel('dV/dt (V/sec)')
        xlabel('Time (ms)')
        xlim([.5 1.75])
        legend('dVm/dt','threshold');
        thisLegend = legend;
        thisLegend.Box = 'off';
        thisAx = gca;
        setAx(thisAx);
    end
    
    %Fill Matrix
    %sweep x spk number; -1s designate no spike
    if isempty(tThresh)~=1
        if length(tThresh)==size(thVals,2) %if noSpks(thisSweep) == maxSpks(1stSweep:thisSweep-1)
            thVals(ii,:)=file(tThresh,1,ii);
            thTimes(ii,:)=tThresh;
        elseif length(tThresh)<size(thVals,2) %if noSpks(thisSweep) < maxSpks(1stSweep:thisSweep-1)
            %adds -1 to complete row; script will ignore -1s
            noRows = size(thVals,2)-length(tThresh);
            addrows=-1.*ones(1,noRows);
            thVals(ii,:)=vertcat(file(tThresh,1,ii),addrows');
            tThresh=horzcat(tThresh,addrows);
            thTimes(ii,:)=tThresh;
        elseif length(tThresh)>size(thVals,2) %if noSpks(thisSweep) > maxSpks(1stSweep:thisSweep-1)
            %adds -1 to complete all other rows; script will ignore -1s
            noCol = length(tThresh)-size(thVals,2);
            addcol=-1.*ones(size(thVals,1),noCol);
            thVals=horzcat(thVals,addcol);
            thTimes=horzcat(thTimes,addcol);
            thVals(ii,:)=file(tThresh,1,ii);
            thTimes(ii,:)=tThresh;
        end
    elseif isempty(tThresh)==1
        thVals(ii,:)=-1;
        thTimes(ii,:)=-1;
    end
    clf
end
close(figure(1))

end

function Rm =getRm(file,blstart,blend,injstart,fs,rheobase,dI,cRange)
%uses all sweeps prior to rheobase, first 100ms of injection compared to 100ms pre injection
if rheobase > 3 %if you start with negative current this shouldn't be a problem
    
    %IV plot to get membrane resistance
    %get dI and dV
    dI(rheobase:end) = [];
    dV=-1.*ones(1,rheobase-1);
    for ii = 1:rheobase-1
        %take voltage from 50-150ms post injection, 50ms for membrane charging & go for 100ms before potential HCN activation (HCN would introduce non-linearities)
        dV(ii) = mean(file(injstart+0.05*fs:(injstart+0.15*fs-1),1,ii))-mean(file(blstart:blend,1,ii));
    end
    
    %get linear fit of IV plot
    IVfit=nlinfit(dI',dV,@dr_fitline,[0 1]);
    plotIVfit=IVfit(1).*dI+IVfit(2);
    Rm=1000*IVfit(1); %gives Rm in MOhms
    
    %Plot IV
    rmFig = figure(2);
    rmFig.Position = [20 475 275 210];
    hold on
    scatter(dI,dV,50,cRange(1,:),'fill')
    plot(dI,plotIVfit,'linewidth',2,'color',cRange(end,:))
%     ylim([-40 20])
    ylabel('Vm (mV)')
    xlabel('current (pA)')
    title({'IV plot; Rm = ' num2str(Rm) ' M\Omega'})
    ivAx = gca;
    setAx(ivAx);
    ivAx.XAxisLocation = 'origin'; ivAx.YAxisLocation = 'origin';
else
    Rm = NaN;
end
end

function [mtau] = getTau(data,sweep,injstart,samplerate,cRange)
%set up single exp fit
xvals=(injstart:injstart+.25*samplerate)'; %first 250ms of injection
endFallingV=find(data(xvals,1,sweep)==min(data(xvals,1,sweep)))+xvals(1); %when voltage reaches most negative Vm
yvals=(data(xvals(1):endFallingV(1),1,sweep));
yPrime = data(xvals,1,sweep);
newXVals=(xvals(1):endFallingV(1))';
dt=newXVals-newXVals(1);
dtPrime = xvals-xvals(1);
dV=yvals-yvals(end);
dVPrime = yPrime - yvals(end);
betaBest=nlinfit(dt,dV,@exp1fit,[dV(1) mean(dt)]);

mtau=betaBest(2)/10; %gets tau_m in milliseconds

%plot
tauFig=figure(3);
tauFig.Position = [305 475 275 210];
hold on
plot(dtPrime,dVPrime,'color',cRange(sweep,:),'linewidth',2)
plot(dt,exp1fit(betaBest,dt),'color',cRange(end,:),'linewidth',2,'linestyle','-.')
title({'Membrane \tau = ' mtau ' ms'})
tauAx = gca;
setAx(tauAx);
tauAx.XTick = 0:.05*samplerate:.25*samplerate;
tauAx.XTickLabel = {'0' '50' '100' '150' '200' '250'};
xlabel('time (ms)')
ylabel('relative Vm (mV)')


end

function vsag = getSag(data,sweep,injstart,injend,blstart,blend,samplerate,cRange)
%returns voltage sag
% 100*(minV-steadystateV)/(minV-baselineV)
% minV is the value of the minumum point in the hyperpolarizing injection
% steadystateV is the mean value over the last 200ms of the hyperpolarizing injection
% baselineV is the mean value over the first 500ms of the sweep

t = 1/samplerate:1/samplerate:size(data,1)/samplerate;
minV=min(data(injstart:injend,1,sweep));
ssV=mean(data((injend-0.2*samplerate):injend,1,sweep));
blV=mean(data(blstart:blend,1,sweep));

vsag=100*((minV-ssV)/(minV-blV));

%plot
sagFig=figure(4);
sagFig.Position = [20 175 275 210];
hold on
plot(t,data(:,1,sweep),'color',cRange(1,:),'linewidth',2)
line([injstart./samplerate-.05 injstart./samplerate-.05],[blV minV],'linewidth',2','color',cRange(round(.67*size(cRange,1)),:))
line([injstart./samplerate-.025 injstart./samplerate-.025],[blV ssV],'linewidth',2','color',cRange(end,:))
text(injstart./samplerate+.025, blV,'maximum voltage response','color',cRange(round(.67*size(cRange,1)),:),'FontWeight','bold')
text(injstart./samplerate+.025, blV-5,'steady state voltage response','color',cRange(end,:),'FontWeight','bold')
xlim([injstart./samplerate-.1 injend./samplerate+.15])
title({'Voltage sag = ' num2str(vsag) '%'})
sagAx = gca;
setAx(sagAx);
xlabel('time (s)')
ylabel('Vm (mV)')
end

function [peaks, times, number] = getRbndspks(data,sweep,t1,stop,cRange,samplerate)
%gets rebound spikes, define as voltage peaks above -15mV

t = 1/samplerate:1/samplerate:size(data,1)/samplerate;

[allpeaks,alltimes]=findpeaks(data(t1:stop,1,sweep));

%plot
rbndFig = figure(5);
rbndFig.Position = [305 175 275 210];
hold on
plot(t,data(:,1,sweep),'color',cRange(sweep,:),'linewidth',2)
scatter(alltimes./samplerate,allpeaks,'markeredgecolor',cRange(end,:))
xlim([t1./samplerate t1./samplerate+1])
setAx(gca);
xlabel('time (s)')
ylabel('Vm (mv)')

peaks=allpeaks(allpeaks>-15);
times=alltimes(allpeaks>-15)-1+t1;
number=length(peaks);
title({num2str(number) ' rebound spikes'})

end

function [maxFiring, maxsweep, isi_cv, isi_sweep, first_isi, second_isi, isi_fr, delay] = getMaxFiring(data,noSweeps,rheoSweep,spkMatrix,samplerate,injStart,cRange)
t = 1/samplerate:1/samplerate:size(data,1)/samplerate;

%find last current injection before AP reduction
nospks=zeros(noSweeps,1);
for ii = rheoSweep:noSweeps
    nospks(ii)=sum(spkMatrix(ii,:)~=-1);
end
sweeps_with_multispikes = find(nospks>3);
isi_sweep = rheoSweep+1;

dspks=zeros(length(rheoSweep:noSweeps),1);
for ii = rheoSweep:noSweeps
    dspks(ii)=nospks(ii)-nospks(ii-1);
end
negsweeps=find(dspks<0);
if isempty(negsweeps) ~= 1
    maxsweep=negsweeps(1)-1;
else
    maxsweep=noSweeps;
end

%get isi during first 200ms of spiking
tFirstSpk = spkMatrix(maxsweep,1);
spks=(spkMatrix(maxsweep,:)>=tFirstSpk)&(spkMatrix(maxsweep,:)<=(tFirstSpk+.2*samplerate));
tspks=spkMatrix(maxsweep,spks); %gets time for spk threshold

%get isi for first multi-spike sweep for isi
tFirstSpk_multi = spkMatrix(isi_sweep,1);
spks_isi=(spkMatrix(isi_sweep,:)>=tFirstSpk_multi);
tspks_isi=spkMatrix(isi_sweep,spks_isi); %gets time for spk threshold

%get ISI in sec
for ii = 1:length(tspks)
    if ii ~=length(tspks)
        ISI(ii)=(tspks(ii+1)-tspks(ii))/samplerate;
    else
        if ii ~= length(find(spkMatrix(maxsweep,:)~=-1))
            ISI(ii)=(spkMatrix(maxsweep,ii+1)-tspks(ii))/samplerate;
        end
    end
end
meanISI=mean(ISI);

%get ISI from first sweep after rheobase in ms
isi_multi = diff(tspks_isi)./samplerate;
isi_multi = 1000.*isi_multi;
isi_cv = std(isi_multi)./mean(isi_multi);
first_isi = isi_multi(1); second_isi = isi_multi(2);
isi_fr = 1/((mean(isi_multi))/1000);

%get delay bias
tspk_sec = (tspks_isi-injStart)./samplerate;
n_spks = length(tspk_sec);
early_spk = sum(tspk_sec<=0.5); late_spk = sum(tspk_sec>0.5);
delay = (late_spk - early_spk)/n_spks;

%get FR in Hz
maxFiring=1/meanISI;

%plot
maxfrFig = figure(6);
maxfrFig.Position = [650 235 620 400];
hold on
plot(t,data(:,1,maxsweep),'linewidth',2,'color',cRange(maxsweep,:))
xlim([injStart./samplerate-.05 injStart./samplerate+1.1])
title({'Max Firing Rate = ' num2str(maxFiring) ' Hz'})
xlabel('time (s)')
ylabel('Vm (mV)')
setAx(gca);

end
function [curve] = getFICurve(data,noSweeps,rheoSweep,sweepDelta,spkMatrix,samplerate,injStart,cRange)
%gets and plots current injection/firing rate relationship
curve = -1.*ones(noSweeps,2);
t = 1/samplerate:1/samplerate:size(data,1)/samplerate;

for ii = 1:noSweeps
    if ii < rheoSweep
        curve(ii,1) = 0;
        curve(ii,2) = sweepDelta(ii,2);
    elseif ii >= rheoSweep
        spks = []; tspks = [];
        %get isi
        spks=(spkMatrix(ii,:)>0);
        tspks=spkMatrix(ii,spks); %gets time for spk threshold
        
        %get isi in sec
        isi = [];
        for jj = 1:length(tspks)
            if jj ~= length(tspks)
                isi(jj) = (tspks(jj+1)-tspks(jj))/samplerate;
            else
                if jj ~= length(spkMatrix(ii,:))
                    if spkMatrix(ii,jj+1) ~= -1
                        isi(jj) = (spkMatrix(ii,jj+1)-tspks(jj))/samplerate;
                    end
                end
            end
        end
        if isempty(isi) ~= 1
            mISI = mean(isi);
            %get FR in Hz
            curve(ii,1) = mISI^(-1);
        else
            curve(ii,1) = 0;
        end
        curve(ii,2) = sweepDelta(ii,2);
    end
end

deleteThese = find(curve(:,2)<0);
curve(deleteThese,:) = [];


%plot
fiFig = figure(8);
fiFig.Position = [1300 225 300 400];
hold on
plot(curve(:,2),curve(:,1),'linewidth',2,'color',cRange(1,:))
title({'FR plot'})
xlabel('current injection (pA)')
ylabel('Sustained Firing Rate (Hz)')
setAx(gca);

end

function [peakval, tpeak, nopeaks] = getPeaks(data,threshVals,threshTimes,sweep)
%detects action potential peaks

nopeaks=sum(threshVals(sweep,:)~=-1);
peakval=zeros(nopeaks,1);
tpeak=zeros(nopeaks,1);
if nopeaks~=0 %if APs in current sweep
    for ii=1:nopeaks
        if ii ~= nopeaks
            [peakval(ii),tpeak(ii)]=max(data(threshTimes(sweep,ii):threshTimes(sweep,ii+1),1,sweep)); %finds max point between threshold_n and threshold_n+1
        elseif ii == nopeaks
            [peakval(ii),tpeak(ii)]=max(data(threshTimes(sweep,ii):threshTimes(sweep,ii)+100,1,sweep)); %finds max point between threshold_n and threshold_n+1
        end
        tpeak(ii)=tpeak(ii)+threshTimes(sweep,ii)-1; %fixes time
    end
end
end

function [AHP, tAHP, ahp_tau] = getAHP(data,peaktimes,threshtimes,injend,sweep,samplerate)
%detects AHP, defined as the largest hyperpolarization deflection within the first 15ms after AP
AHP=zeros(length(peaktimes),1);
tAHP=zeros(length(peaktimes),1);
for ii = 1:length(peaktimes)
    if ii ~=length(peaktimes)
        deltaspk=threshtimes(ii+1)-threshtimes(ii);
        if deltaspk <= 0.015*samplerate
            [AHP(ii), tAHP(ii)]=min(data(peaktimes(ii):threshtimes(ii+1),1,sweep));
        else
            [AHP(ii), tAHP(ii)]=min(data(peaktimes(ii):threshtimes(ii)+0.015*samplerate,1,sweep));
        end
        tAHP(ii)=tAHP(ii)-1+peaktimes(ii); %corrects AHP time
        AHP_trace(ii,:) = data(tAHP(ii):tAHP(ii)+0.05*samplerate,1,sweep);
    else
        deltaend=injend-threshtimes(ii);
        if deltaend <= 0.1*samplerate
            if deltaend <=0
                AHP(ii)=[];
                tAHP(ii)=[];
            else
                if (injend-peaktimes(ii)) <= 0.02*samplerate
                    AHP(ii)=[];
                    tAHP(ii)=[];
                else
                    [AHP(ii), tAHP(ii)]=min(data(peaktimes(ii):injend,1,sweep));
                    tAHP(ii)=tAHP(ii)-1+peaktimes(ii); %corrects AHP time
                end
            end
        else
            [AHP(ii), tAHP(ii)]=min(data(peaktimes(ii):threshtimes(ii)+0.015*samplerate,1,sweep));
            tAHP(ii)=tAHP(ii)-1+peaktimes(ii); %corrects AHP time
        end
    end
end

%set up single exp fit for AHP decay
mean_trace = mean(AHP_trace);
xvals=(1:.03*samplerate)'; %first 10-50ms after peak
endRisingV=find(mean_trace(xvals)==max(mean_trace(xvals)))+xvals(1); %when voltage reaches most pos Vm
yvals=(mean_trace(xvals(1):endRisingV(1)));
yPrime = mean_trace(xvals);
newXVals=(xvals(1):endRisingV(1))';
dt=newXVals-newXVals(1);
dtPrime = xvals-xvals(1);
dV=yvals-yvals(end);
dVPrime = yPrime - yvals(end);
betaBest=nlinfit(dt,dV',@exp1fit,[dV(1) mean(dt)]);

ahp_tau=betaBest(2)/10; %gets tau_m in milliseconds

%plot
tauFig=figure(30);
tauFig.Position = [305 475 275 210];
hold on
plot(dtPrime,dVPrime,'linewidth',2)
plot(dt,exp1fit(betaBest,dt),'linewidth',2,'linestyle','-.')
title({'AHP \tau = ' ahp_tau ' ms'})
tauAx = gca;
setAx(tauAx);
tauAx.XTick = 0:.01*samplerate:.05*samplerate;
tauAx.XTickLabel = {'0' '10' '20' '30' '40' '50'};
xlabel('time (ms)')
ylabel('relative Vm (mV)')

end

function [hw,halfamp1,tHalf1,halfamp2,tHalf2] = getHalfwidth(data,sweep,samplerate,threshtimes,peaktime,amp)
%defined as the width (in msec) of the AP waveform at the half-amplitude
halfamp1=zeros(length(threshtimes),1);
halfamp2=zeros(length(threshtimes),1);
hw=zeros(length(threshtimes),1);
tHalf1=zeros(length(threshtimes),1);
tHalf2=zeros(length(threshtimes),1);
for ii = 1:length(threshtimes)
    dfrom_halfamp1=[];
    dfrom_halfamp2=[]; %#ok<*NASGU>
    halfamp=amp(ii)/2;
    dfrom_halfamp1=abs(data(threshtimes(ii):peaktime(ii)-1,1,sweep)-(data(peaktime(ii),1,sweep)-halfamp));
    [~,t_mindist1]=min(dfrom_halfamp1);
    dfrom_halfamp2=abs(data(peaktime(ii):peaktime(ii)+0.0025*samplerate,1,sweep)-(data(peaktime(ii),1,sweep)-halfamp));
    [~,t_mindist2]=min(dfrom_halfamp2);
    tHalf1(ii)=threshtimes(ii)+t_mindist1-1;
    halfamp1(ii)=data(tHalf1(ii),1,sweep);
    tHalf2(ii)=peaktime(ii)+t_mindist2-1;
    halfamp2(ii)=data(tHalf2(ii),1,sweep);
    hw(ii)=(tHalf2(ii)-tHalf1(ii))/(samplerate/1000);
end
end

function [ratio,std_isi] = getRatioFR(threshtimes,sweep,nospks)
%ratio of last 3 ISI to first ISI in sweep

tspks=threshtimes(sweep,1:nospks-1); %gets time for spk threshold
%get ISI in dps
for ii = 1:length(tspks)
    if ii ~=length(tspks)
        ISI(ii)=(tspks(ii+1)-tspks(ii));
    else
        ISI(ii)=(threshtimes(sweep,ii+1)-tspks(ii));
    end
end

if length(ISI) >= 3
    ratio=ISI(1)/mean([ISI(end) ISI(end-1) ISI(end-2)]);
    std_isi = std(ISI./10); %gives std in ms
else
    ratio=NaN;
end
end