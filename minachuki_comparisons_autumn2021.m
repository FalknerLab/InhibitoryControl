%created 10-13-21 %modified 10-13-21 %author eartha mae guthman

%edit analysis of mean responses to mean sweep rather than aligned by max rise

close all
clearvars

%% INITS
inputLabels = {'mpoa','shell'};

%colors
cc.mpoa = [180 151 214]./255; %wisteria
cc.shell = .8.*[89 217 153]./255; %medium aquamarine
colormap(brewermap(256,'Spectral'))
spectrals = colormap(brewermap(256,'Spectral'));

%single stim
noCell.oIPSC.mpoa = 13;
noCell.oIPSC.shell = 11;
rowID.oIPSC.mpoa = 34;
rowID.oIPSC.shell = 26;
singleGrpStr = {'mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','shell','shell','shell','shell','shell','shell','shell','shell','shell','shell','shell'};

%paired pulse
noCell.pp.mpoa = 11;
noCell.pp.shell = 12;
rowID.pp.mpoa = 45;
rowID.pp.shell = 42;
ppGrpStr = {'mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','mpoa','shell','shell','shell','shell','shell','shell','shell','shell','shell','shell', 'shell', 'shell'};

%membrane props
noCell.mem_props = 51;
grouping = []; %group IDs based off of hierarchical clustering/ward's method
grouping(1:35) = 1; grouping(36:51) = 2;


%% LOAD DATA
%google sheet
sheetID = '1s3Wbv1VEYNm76z7HRIZSV8yjnmek1ZWt62ig8qeOiVg'; %google sheet ID for ephys data
sheet_name = 'My Sheet';
url_name = sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',sheetID, sheet_name);
sheet_data = webread(url_name);

%single stim
data.oIPSC.mpoa = table2array(sheet_data(rowID.oIPSC.mpoa,2:noCell.oIPSC.mpoa+1));
data.oIPSC.mpoa_ap = table2array(sheet_data(rowID.oIPSC.mpoa+6,2:noCell.oIPSC.mpoa+1));
data.oIPSC.shell = table2array(sheet_data(rowID.oIPSC.shell,2:noCell.oIPSC.shell+1));
data.oIPSC.shell_ap = table2array(sheet_data(rowID.oIPSC.shell+6,2:noCell.oIPSC.shell+1));

%paired pulse
data.pp.mpoa = table2array(sheet_data(rowID.pp.mpoa,2:noCell.pp.mpoa+1));
data.pp.mpoa_ap = table2array(sheet_data(rowID.pp.mpoa+1,2:noCell.pp.mpoa+1));
data.pp.shell = table2array(sheet_data(rowID.pp.shell,2:noCell.pp.shell+1));
data.pp.shell_ap = table2array(sheet_data(rowID.pp.shell+1,2:noCell.pp.shell+1));

%membrane properties
data.mem_props = table2array(sheet_data(2:23,2:noCell.mem_props+1));
labels.mem_props = table2cell(sheet_data(2:23,1));

%% ANALYZE DATA
%single stim
%normality tests
for ii = 1:length(inputLabels)
    [holdH.oIPSC.(inputLabels{ii}),~] = swtest(data.oIPSC.(inputLabels{ii}),0.05,1);
end
h_oIPSC = sum([holdH.oIPSC.(inputLabels{1}) holdH.oIPSC.(inputLabels{2})]);

%comparison tests
if h_oIPSC == 0 %if data are normal, unpaired ttest
    [h_comparisons.oIPSC.unpaired_ttest, p.oIPSC.unpaired_ttest] = ttest2(data.oIPSC.mpoa,data.oIPSC.shell);
    p.oIPSC.bartlett = vartestn([data.oIPSC.mpoa'; data.oIPSC.shell'],...
        singleGrpStr,'TestType','Bartlett'); %bartlett test for homogeneity of variances
else %if data are not normal, mann whitney u test
    [p.oIPSC.mwu_test, h_comparisons.oIPSC.mwu_test] = ranksum(data.oIPSC.mpoa,data.oIPSC.shell);
    p.oIPSC.brown_forsythe = vartestn([data.oIPSC.mpoa'; data.oIPSC.shell'],...
        singleGrpStr,'TestType','BrownForsythe'); %brown forsythe test for homogeneity of variances
end

%summary stats
if h_oIPSC == 0
    for ii = 1:length(inputLabels)
        means.oIPSC.(inputLabels{ii}) = mean(data.oIPSC.(inputLabels{ii}));
        sem_error.oIPSC.(inputLabels{ii}) = sem(data.oIPSC.(inputLabels{ii}),2);
    end
else
    for ii = 1:length(inputLabels)
        medians.oIPSC.(inputLabels{ii}) = median(data.oIPSC.(inputLabels{ii}));
        [lower_error.oIPSC.(inputLabels{ii}), upper_error.oIPSC.(inputLabels{ii})] = iqrError(data.oIPSC.(inputLabels{ii}),2);
    end
end

for ii = 1:length(inputLabels)
    [holdH.pp.(inputLabels{ii}),~] = swtest(data.pp.(inputLabels{ii}),0.05,1);
end
h_pp = sum([holdH.pp.(inputLabels{1}) holdH.pp.(inputLabels{2})]);

%relationship to AP position?
lr.oIPSC.shell = fitlm(data.oIPSC.shell_ap,data.oIPSC.shell);
lr.oIPSC.mpoa = fitlm(data.oIPSC.mpoa_ap,data.oIPSC.mpoa);

%comparison tests
if h_pp == 0 %if data are normal, unpaired ttest
    [h_comparisons.pp.unpaired_ttest, p.pp.unpaired_ttest] = ttest2(data.pp.mpoa,data.pp.shell);
    p.pp.bartlett = vartestn([data.pp.mpoa'; data.pp.shell'],...
        ppGrpStr,'TestType','Bartlett'); %bartlett test
else %if data are not normal, mann whitney u test
    [p.pp.mwu_test, h_comparisons.pp.mwu_test] = ranksum(data.pp.mpoa,data.pp.shell);
    p.pp.brown_forsythe = vartestn([data.pp.mpoa'; data.pp.shell'],...
        ppGrpStr,'TestType','BrownForsythe'); %brown forsythe test for homogeneity of variances
end


%summary stats
if h_pp == 0
    for ii = 1:length(inputLabels)
        means.pp.(inputLabels{ii}) = mean(data.pp.(inputLabels{ii}));
        sem_error.pp.(inputLabels{ii}) = sem(data.pp.(inputLabels{ii}),2);
    end
else
    for ii = 1:length(inputLabels)
        medians.pp.(inputLabels{ii}) = median(data.pp.(inputLabels{ii}));
        [lower_error.pp.(inputLabels{ii}), upper_error.pp.(inputLabels{ii})] = iqrError(data.pp.(inputLabels{ii}),2);
    end
end

%relationship to AP position?
lr.pp.shell = fitlm(data.pp.shell_ap,data.pp.shell);
lr.pp.mpoa = fitlm(data.pp.mpoa_ap,data.pp.mpoa);

%membrane properties
%normalize values
data.norm_props = -1.*ones(size(data.mem_props));
for ii = 1:size(data.mem_props,1)
    data.norm_props(ii,:) = zscore(data.mem_props(ii,:));
end

%clustering
%PCA to hierarchical
[data.pca_coeff,data.pca_score,data.pca_latent] = pca(data.norm_props');
var_thresh = 90; %variance to select PCA for clustering

%get PCA latents as % var
data.pca_var = 100.*(data.pca_latent/sum(data.pca_latent));
data.cum_var = cumsum(data.pca_var);
cutoff = find(data.cum_var>var_thresh);
try
    cutoff(2:end) = [];
end
data.pca_matrix = data.pca_score(:,1:cutoff);

%hierarchical cluster
Z.pca = linkage(data.pca_matrix,'ward');
[~,~,outperm] = dendrogram(Z.pca,0);
groups = [];
groups(outperm) = grouping;

%tSNE
data.mem_props_tSNE = tsne(data.pca_matrix,'Perplexity',8,'Distance','cosine');

%normality tests
h_mem_props = zeros(1,size(data.mem_props,1));
for ii = 1:size(data.mem_props,1)
    for clus = 1:2
        [holdH.mem_props(clus,ii),~] = swtest(data.mem_props(ii,groups==clus),0.05,1);
    end
    h_mem_props(ii) = sum([holdH.mem_props(1,ii) holdH.mem_props(2,ii)]);

    %comparison tests
    if h_mem_props(ii) == 0 %if data are normal, unpaired ttest
        [h_comparisons.mem_props(ii).center, p.mem_props(ii).center, ~, stats.mem_props(ii).center] = ttest2(data.mem_props(ii,groups==1),data.mem_props(ii,groups==2));
        p.mem_props(ii).var = vartestn([data.mem_props(ii,groups==1)';data.mem_props(ii,groups==2)'],...
            grouping,'TestType','Bartlett'); %bartlett test for homogeneity of variances
    else %if data are not normal, mann whitney u test
        [p.mem_props(ii).center, h_comparisons.mem_props(ii).center, stats.mem_props(ii).center] = ranksum(data.mem_props(ii,groups==1),data.mem_props(ii,groups==2));
        p.mem_props(ii).var = vartestn([data.mem_props(ii,groups==1)';data.mem_props(ii,groups==2)'],...
            grouping,'TestType','BrownForsythe'); %brown forsythe test for homogeneity of variances
    end

    %summary stats
    if h_mem_props(ii) == 0
        for clus = 1:2
            means.mem_props(ii,clus) = mean(data.mem_props(ii,groups==clus));
            sem_error.mem_props(ii,clus) = sem(data.mem_props(ii,groups==clus),2);
        end
    else
        for clus = 1:2
            medians.mem_props(ii,clus) = median(data.mem_props(ii,groups==clus));
            [lower_error.mem_props(ii,clus), upper_error.mem_props(ii,clus)] = iqrError(data.mem_props(ii,groups==clus),2);
        end
    end
end

%% PLOT DATA
close all

%single stim
%fig and labels
oStimFig = figure;
oStimFig.Position = [-2500 750 340 420];
hold on
title('oIPSC stimulation')
ylabel('oIPSC amplitude')
xlim([0.5 2.5])
ylim([1.1*min([min(data.oIPSC.mpoa) min(data.oIPSC.shell)]) 0])
oStimAx = gca;
setAx(oStimAx);
oStimAx.FontSize = 24;

%scatter
scatter_var = 0.05;
for ii = 1:length(inputLabels)
    scatter_pos.oIPSC.(inputLabels{ii}) = scatter_var.*randn(noCell.oIPSC.(inputLabels{ii}),1) + ii;
    h.scatter.oIPSC = scatter(scatter_pos.oIPSC.(inputLabels{ii}),data.oIPSC.(inputLabels{ii}),150,cc.(inputLabels{ii}),'filled');
    if h_oIPSC == 0
        errorbar(ii+.2,means.oIPSC.(inputLabels{ii}),sem_error.oIPSC.(inputLabels{ii}),...
            'marker','o','markersize',20,'markerfacecolor',cc.(inputLabels{ii}),'MarkerEdgeColor','none',...
            'CapSize',0,'linewidth',6,'Color',cc.(inputLabels{ii}))
    else
        errorbar(ii+.33,medians.oIPSC.(inputLabels{ii}),lower_error.oIPSC.(inputLabels{ii}),upper_error.oIPSC.(inputLabels{ii}),...
            'marker','o','markersize',20,'markerfacecolor',cc.(inputLabels{ii}),'MarkerEdgeColor','none',...
            'CapSize',0,'linewidth',6,'Color',cc.(inputLabels{ii}))
    end
end
oStimAx.XTickLabels = {['n = ' num2str(noCell.oIPSC.(inputLabels{1}))],['n = ' num2str(noCell.oIPSC.(inputLabels{2}))]};

%additional aesthetics
oStimLgd = legend(inputLabels{1},'',inputLabels{2},'','Box','off','Location','best','FontSize',20);

%linear regression
oIPSC_lrFig_shell = figure;
oIPSC_lrFig_shell.Position = [-1900 750 240 420];
hold on
lr_plot_oIPSC.shell = plot(lr.oIPSC.shell);
lr_plot_oIPSC.shell(1).Marker = 'o'; 
lr_plot_oIPSC.shell(1).MarkerFaceColor = cc.shell; lr_plot_oIPSC.shell(1).MarkerEdgeColor = 'none';
lr_plot_oIPSC.shell(2).Color = cc.shell; lr_plot_oIPSC.shell(2).LineWidth = 2;
lr_plot_oIPSC.shell(3).Color = cc.shell; lr_plot_oIPSC.shell(3).LineWidth = 1;
lr_plot_oIPSC.shell(4).Color = cc.shell; lr_plot_oIPSC.shell(4).LineWidth = 1;
title('single stim')
ylabel('oIPSC amplitude (pA)')
xlabel('AP position (mm)')
ylim([1.1*min([min(data.oIPSC.mpoa) min(data.oIPSC.shell)]) 0])
xlim([-2 -1.3])
lr_shellAx = gca;
setAx(lr_shellAx); lr_shellAx.Legend.Visible = 'off';
lr_shellAx.FontSize = 24;
lr_shellAx.XDir = 'reverse';

oIPSC_lrFig_mpoa = figure;
oIPSC_lrFig_mpoa.Position = [-2140 750 240 420];
hold on
lr_plot_oIPSC.mpoa = plot(lr.oIPSC.mpoa);
lr_plot_oIPSC.mpoa(1).Marker = 'o'; 
lr_plot_oIPSC.mpoa(1).MarkerFaceColor = cc.mpoa; lr_plot_oIPSC.mpoa(1).MarkerEdgeColor = 'none';
lr_plot_oIPSC.mpoa(2).Color = cc.mpoa; lr_plot_oIPSC.mpoa(2).LineWidth = 2;
lr_plot_oIPSC.mpoa(3).Color = cc.mpoa; lr_plot_oIPSC.mpoa(3).LineWidth = 1;
lr_plot_oIPSC.mpoa(4).Color = cc.mpoa; lr_plot_oIPSC.mpoa(4).LineWidth = 1;
title('single stim')
ylabel('oIPSC amplitude (pA)')
xlabel('AP position (mm)')
ylim([1.1*min([min(data.oIPSC.mpoa) min(data.oIPSC.shell)]) 0])
xlim([-2 -1.3])
lr_mpoaAx = gca; lr_mpoaAx.Legend.Visible = 'off';
setAx(lr_mpoaAx);
lr_mpoaAx.FontSize = 24;
lr_mpoaAx.XDir = 'reverse';

%paired pulse
%fig and labels
ppFig = figure;
ppFig.Position = [-2500 1 340 420];
hold on
title('paired pulse')
ylabel('pp ratio (peak_2 / peak_1)')
xlim([0.5 2.5])
ylim([0 1.1*max([max(data.pp.mpoa) max(data.pp.shell)])])
ppAx = gca;
setAx(ppAx);
ppAx.FontSize = 24;

%scatter
scatter_var = 0.05;
for ii = 1:length(inputLabels)
    scatter_pos.pp.(inputLabels{ii}) = scatter_var.*randn(noCell.pp.(inputLabels{ii}),1) + ii;
    h.scatter.pp = scatter(scatter_pos.pp.(inputLabels{ii}),data.pp.(inputLabels{ii}),150,cc.(inputLabels{ii}),'filled');
    if h_pp == 0
        errorbar(ii+.2,means.pp.(inputLabels{ii}),sem_error.pp.(inputLabels{ii}),...
            'marker','o','markersize',20,'markerfacecolor',cc.(inputLabels{ii}),'MarkerEdgeColor','none',...
            'CapSize',0,'linewidth',6,'Color',cc.(inputLabels{ii}))
    else
        errorbar(ii+.33,medians.pp.(inputLabels{ii}),lower_error.pp.(inputLabels{ii}),upper_error.pp.(inputLabels{ii}),...
            'marker','o','markersize',20,'markerfacecolor',cc.(inputLabels{ii}),'MarkerEdgeColor','none',...
            'CapSize',0,'linewidth',6,'Color',cc.(inputLabels{ii}))
    end
end
ppAx.XTickLabels = {['n = ' num2str(noCell.pp.(inputLabels{1}))],['n = ' num2str(noCell.pp.(inputLabels{2}))]};

%additional aesthetics
ppLgd = legend(inputLabels{1},'',inputLabels{2},'','Box','off','Location','best','FontSize',20);

%linear regression
pp_lrFig_shell = figure;
pp_lrFig_shell.Position = [-1900 1 240 420];
hold on
lr_plot_pp.shell = plot(lr.pp.shell);
lr_plot_pp.shell(1).Marker = 'o'; 
lr_plot_pp.shell(1).MarkerFaceColor = cc.shell; lr_plot_pp.shell(1).MarkerEdgeColor = 'none';
lr_plot_pp.shell(2).Color = cc.shell; lr_plot_pp.shell(2).LineWidth = 2;
lr_plot_pp.shell(3).Color = cc.shell; lr_plot_pp.shell(3).LineWidth = 1;
lr_plot_pp.shell(4).Color = cc.shell; lr_plot_pp.shell(4).LineWidth = 1;
title('paired pulse')
ylabel('pp ratio (peak_2 / peak_1)')
xlabel('AP position (mm)')
ylim([0 1.1*max([max(data.pp.mpoa) max(data.pp.shell)])])
xlim([-2 -1.3])
lr_shellAx = gca;
setAx(lr_shellAx); lr_shellAx.Legend.Visible = 'off';
lr_shellAx.FontSize = 24;
lr_shellAx.XDir = 'reverse';

pp_lrFig_mpoa = figure;
pp_lrFig_mpoa.Position = [-2140 1 240 420];
hold on
lr_plot_pp.mpoa = plot(lr.pp.mpoa);
lr_plot_pp.mpoa(1).Marker = 'o'; 
lr_plot_pp.mpoa(1).MarkerFaceColor = cc.mpoa; lr_plot_pp.mpoa(1).MarkerEdgeColor = 'none';
lr_plot_pp.mpoa(2).Color = cc.mpoa; lr_plot_pp.mpoa(2).LineWidth = 2;
lr_plot_pp.mpoa(3).Color = cc.mpoa; lr_plot_pp.mpoa(3).LineWidth = 1;
lr_plot_pp.mpoa(4).Color = cc.mpoa; lr_plot_pp.mpoa(4).LineWidth = 1;
title('paired pulse')
ylabel('pp ratio (peak_2 / peak_1)')
xlabel('AP position (mm)')
ylim([0 1.1*max([max(data.pp.mpoa) max(data.pp.shell)])])
xlim([-2 -1.3])
lr_mpoaAx = gca; lr_mpoaAx.Legend.Visible = 'off';
setAx(lr_mpoaAx);
lr_mpoaAx.FontSize = 24;
lr_mpoaAx.XDir = 'reverse';

%membrane properties
mem_propsFig = figure;
mem_propsFig.Position = [1 1 1100 1016];
hold on
subplot(4,6,1:2)
h.dendro = dendrogram(Z.pca,0);
dendroAx = gca;
setAx(dendroAx);
dendroAx.YLim(1) = 0; dendroAx.XTick = [];

scatter_var = 0.1; 
for ii = 1:size(data.norm_props,1) %loops thru membrane variables
    scatter_pos.memprops = scatter_var.*randn(noCell.mem_props,1) + groups';
    %set up subplot
    subplot(4,6,ii+2)
    gscatter(scatter_pos.memprops,data.mem_props(ii,:),groups)
    xlim([.5 2.5])

    %subplot aesthetics
    mem_propsAx(ii) = gca;
    setAx(mem_propsAx(ii));
    mem_propsAx(ii).FontSize = 12; 
    mem_propsAx(ii).Box = 'off';
    mem_propsAx(ii).Legend.Visible = 'off';
    title(labels.mem_props{ii})
end

% tSNE plot of membrane properties
tSNE_fig = figure;
tSNE_fig.Position  = [1375 1 1740 1336];
subplot(4,6,1)
gscatter(data.mem_props_tSNE(:,1),data.mem_props_tSNE(:,2),groups,[],[],30)
tSNE_ax(1) = gca;
setAx(tSNE_ax(1));
tSNE_ax(1).FontSize = 18;
tSNE_ax(1).Box = 'off';


for ii = 1:size(data.norm_props,1) %loops thru membrane variables
    %set up subplot
    subplot(4,6,ii+1)
    scatter(data.mem_props_tSNE(:,1),data.mem_props_tSNE(:,2),50,data.norm_props(ii,:)','filled')
    colormap('viridis')
    colorbar

    %subplot aesthetics
    tSNE_ax(1+ii) = gca;
    setAx(tSNE_ax(1+ii));
    tSNE_ax(ii+1).FontSize = 18;
    tSNE_ax(ii+1).Box = 'off';
    title(labels.mem_props{ii})
end
