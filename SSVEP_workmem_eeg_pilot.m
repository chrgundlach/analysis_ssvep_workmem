%% script to analyze SSVEP pilot data


%% General Definitions
clearvars
p.path_cue=         '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2025_SSVEP_workmem\eeg\epoch\';
p.path_erp=         '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2025_SSVEP_workmem\eeg\epoch_erp\';
p.path_behavior=    '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2025_SSVEP_workmem\behavior\';
p.files =           {'VP01';'VP02';'VP03'};
p.files2use =       [1 2 3];
p.chanlocs_path=    ['C:\Users\psy05cvd\Dropbox\work\matlab\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_1020.epf'];

p.resample=         256;

p.con1idx =             [1 2 3 4 5 6];
p.con1label =           {'cue_left_2t2d'; 'cue_right_2t2d'; 'cue_left_4t'; 'cue_right_4t'; 'cue_left_2t'; 'cue_right_2t'};
% [RDK1 RDK3; RDK4 RDK6] [RDK1 RDK3; RDK4 RDK6] [RDK1 RDK2; RDK4 RDK5] [RDK1 RDK2; RDK4 RDK5] [RDK1; RDK4] [RDK1; RDK4]
p.con2idx =             {[1 3 5];[2 4 6]};
p.con2label =           {'cue_left';'cue_right'};
p.con3idx =             {[1 2];[3 4]; [5 6]};
p.con3label =           {'2t_2d';'4t';'2t'};
p.con_trigger_cue =     [10 20 30 40 50 60];
p.con_trigger_retent =  [11 21 31 41 51 61];

p.cue_epoch =           [-1 1.5];
p.cue_win2an =          {[-1 0];[0.5 1.5]; [-1 1.5]};

p.retent_epoch =        [-0.5 1];

p.FFT_freqres =         2^12;

p.CSD_flag =            1;

p.erp_filter =          [0 14];

%% read in bdfs and behavior
for i_sub = 1:numel(p.files2use)
    %% behavior
    % load behavioral data
    t.files = dir(fullfile(p.path_behavior,sprintf('%s_timing*.mat',p.files{p.files2use(i_sub)}(1:4))));
    [t.val t.idx ]=max([t.files.datenum]);
    behavior = load(fullfile(p.path_behavior,t.files(t.idx).name));

    RDK_all{i_sub}=behavior.RDK.RDK;
    
    %% EEG
    % read in preprocessed files
    EEG_cue = pop_loadset('filename',[p.files{p.files2use(i_sub)} '_e.set'], 'filepath', p.path_cue);
    EEG_erp = pop_loadset('filename',[p.files{p.files2use(i_sub)} '_e.set'], 'filepath', p.path_erp);
   
    % pop_eegplot(EEG_cue,1,1,1)

    %% do csd transform
    if p.CSD_flag==1
        if  i_sub == 1 % calculate CSD matrix
            try
                CSD.chanmat=ExtractMontage('C:\Users\psy05cvd\Dropbox\work\matlab\software\toolboxes\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG_cue.chanlocs.labels}');
            catch
                CSD.chanmat=ExtractMontage('C:\Users\EEG\Documents\MATLAB\christopher\general_functions\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG_cue.chanlocs.labels}');
            end
            [CSD.G,CSD.H] = GetGH(CSD.chanmat);
        end
        fprintf(1,'\n###\ncalculating CSD transform\n###\n')
        for i_tr = 1:EEG_cue.trials
            % csd of raw data
            EEG_cue.data(:,:,i_tr)= CSDTransform(EEG_cue.data(:,:,i_tr), CSD.G, CSD.H);
        end
        for i_tr = 1:EEG_erp.trials
            % csd of raw data
            EEG_erp.data(:,:,i_tr)= CSDTransform(EEG_erp.data(:,:,i_tr), CSD.G, CSD.H);
        end
    end
    
    %% FFT SSVEP analysis
    t.trialindex = arrayfun(@(x) [EEG_cue.event([ismember([EEG_cue.event.type],x)]).epoch], p.con_trigger_cue,'UniformOutput',false);

    % preallocate memory
    if i_sub == 1
        FFT_data_evo=nan(p.FFT_freqres, EEG_cue.nbchan,numel(p.con_trigger_cue),numel(p.cue_win2an));
    end

    %% do windowed FFT analysis
    for i_win = 1:numel(p.cue_win2an)
        %% windowed analysis
        % select data
        EEGt = pop_select(EEG_cue,'time',p.cue_win2an{i_win});
        % surrogate mean data for evoked spectra
        EEGt_m = pop_select(EEGt,'trial',1);
        % detrend
        EEGt = eegF_Detrend(EEGt,[]);
        %pop_eegplot(EEGt,1,1,1)
        
        % loop across conditions
        for i_con = 1:numel(p.con_trigger_cue)
            % evoked
            EEGt_m.data = detrend(mean(EEGt.data(:,:,t.trialindex{i_con}),3)')';
%             EEGt_m.data = mean(EEGt.data(:,:,t.trialindex{i_con}),3);
            FFT_data_evo(:,:,i_con,i_win,i_sub) = squeeze(abs(fft(EEGt_m.data,p.FFT_freqres,2))*2/size(EEGt.data,2))';
        end
    end
    FFT_freqs = ((0:size(FFT_data_evo,1)-1)/size(FFT_data_evo,1)) * EEGt.srate;


    %% epoch for ERP analysis
    % select data
    if ~isempty(p.erp_filter)
        EEG_erp = pop_eegfiltnew(EEG_erp, p.erp_filter(1), p.erp_filter(2), 8*EEG_erp.srate, 0, [], 0);
    end

    EEGt = pop_rmbase(EEG_erp,[-100 0]);

    t.trialindex = arrayfun(@(x) [EEGt.event([ismember([EEGt.event.type],x)]).epoch], p.con_trigger_retent,'UniformOutput',false);

    % preallocate memory
    if i_sub == 1
        ERP_data=nan(EEG_cue.nbchan,EEGt.pnts ,numel(p.con_trigger_cue),numel(numel(p.files2use)));
    end

    % loop across conditions
    for i_con = 1:numel(p.con_trigger_cue)
        ERP_data(:,:,i_con,i_sub) = mean(EEGt.data(:,:,t.trialindex{i_con}),3);
    end
    ERP_times = EEGt.times;
end


%% plot only electrode positions
figure; topoplot([],EEG_cue.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG_cue.chaninfo);

%% do some plotting | FFT spectra
pl.elec2plot = {'POz';'O1';'Oz';'O2';'Iz'};
pl.elec2plot = {'Oz';'O2';'Iz';'I2'}; %!
pl.elec2plot = {'O2';'PO8';'I2'};
pl.elec2plot = {'Iz';'I2'};
pl.elec2plot = {'P7';'PO7';'O1';'Oz';'O2';'PO8';'P8';'P9';'I1';'Iz';'I2';'P10'};
pl.elec2plot = {'O1';'Oz';'O2';'I1';'Iz';'I2'};

% pl.elec2plot = {'P7';'PO7';'Iz';'I2'};         % right stims two cluster
% pl.elec2plot = {'O2';'PO8';'PO4';'I2'};         % left one cluster
% pl.elec2plot = {'PO4';'PO8';'O2'};         % left stims one cluster !
% pl.elec2plot = {'P7';'PO7'};         % right stims two cluster
% pl.elec2plot = {'P7';'PO7';'P9';'O1';'I1';'Iz'};         % right stims two cluster  !

pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({EEG_cue.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.con2plot = [1 2 3 4 5 6]; % [RDK_all{1}.freq]([1 4]) %always present
pl.con2plot = [1 2 3 4]; % [RDK_all{1}.freq]([1 2 4 5]) %always present

pl.sub2plot = 1:numel(p.files2use);
% pl.sub2plot = 2;
pl.time2plot = [1,2];
pl.time2plot = [3];

pl.data = squeeze(mean(FFT_data_evo(:,pl.elec2plot_i,pl.con2plot,pl.time2plot,pl.sub2plot),[2,3,4]));
pl.freq2plot = [RDK_all{1}.freq];

figure;
set(gcf,'Position',[100 100 600 300],'PaperPositionMode','auto')
plot(FFT_freqs,pl.data,'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on;
plot(FFT_freqs,mean(pl.data,2),'Color','k','LineWidth',2)

xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV/cm²')
title(sprintf('Grand Mean evoked spectrum | t = [%1.0f %1.0f]ms | N = %1.0f', ...
    minmax([p.cue_win2an{pl.time2plot}]*1000), numel(pl.sub2plot)),'Interpreter','none')
vline(pl.freq2plot([1, 2, 4, 5]),'k:',{'left','left','right','right'})
% draw topography with electrode positions
h.a1 = axes('position',[0.68 0.68 0.25 0.25],'Visible','off');
topoplot(find(pl.elec2plot_i),EEG_cue.chanlocs,'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',4,1});

%% do some plotting | topographies
pl.freqs = [RDK_all{1}.freq];
pl.freq2plot = 5;
pl.freqpos = {'left';'left';'left';'right';'right';'right'};

pl.sub2plot = 1:numel(p.files2use);
% pl.sub2plot = 1;
pl.time2plot = [1,2];
pl.time2plot = [3];

pl.con2plot = [1 2 3 4 5 6]; % [RDK_all{1}.freq]([1 4]) %always present
pl.con2plot = [1 2 3 4]; % [RDK_all{1}.freq]([1 2 4 5]) %always present

pl.freqrange = [-0.1 0.1];
pl.fidx = dsearchn(FFT_freqs', (pl.freqs(pl.freq2plot)+pl.freqrange)');

pl.data = squeeze(mean(squeeze(FFT_data_evo(pl.fidx(1):pl.fidx(2),:,pl.con2plot,pl.time2plot,pl.sub2plot)),[1,3,4]));
pl.mdata = mean(pl.data,3,'omitnan');

figure;
set(gcf,'Position',[100 100 300 300],'PaperPositionMode','auto')

pl.clim = [0 max(pl.mdata)];
topoplot(pl.mdata, EEG_cue.chanlocs, ...
    'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
    'whitebk','on');
title(sprintf('evoked SSVEP | %s stim\n[%1.0f %1.0f]ms | %1.1f Hz', ...
        pl.freqpos{pl.freq2plot}, minmax([p.cue_win2an{pl.time2plot}]*1000), ...
        pl.freqs(pl.freq2plot)), ...
        "Interpreter","none")
colorbar


%% do some plotting | CDA topo | diff of conditions
pl.time2plot = [200 1000];
pl.time2plot = [500 1000];
pl.time2plot = [1000 1500];
pl.time2plot = [-1000 0];
pl.tidx = dsearchn(ERP_times', pl.time2plot');

pl.sub2plot = 1:numel(p.files2use);
% pl.sub2plot = 1:2;


pl.contrast = p.con3idx([3,1,2]);
pl.contrast_label = p.con3label([3,1,2]);

pl.data = nan([size(ERP_data,1),numel(pl.tidx(1):pl.tidx(2)), size(pl.contrast,1) numel(pl.sub2plot)]);

for i_contr = 1:size(pl.contrast,1)
    pl.data(:,:,i_contr,:)=ERP_data(:,pl.tidx(1):pl.tidx(2),pl.contrast{i_contr}(1),pl.sub2plot) - ERP_data(:,pl.tidx(1):pl.tidx(2),pl.contrast{i_contr}(2),pl.sub2plot);
%     pl.data(:,:,i_contr,:)=ERP_data(:,pl.tidx(1):pl.tidx(2),pl.contrast{i_contr}(1),pl.sub2plot);
%     pl.data(:,:,i_contr,:)=ERP_data(:,pl.tidx(1):pl.tidx(2),pl.contrast{i_contr}(2),pl.sub2plot);
end

pl.mdata = squeeze(mean(pl.data,[2,4]));

figure;
set(gcf,'Position',[100 100 700 180],'PaperPositionMode','auto')
for i_contr = 1:size(pl.mdata,2)
    h.s(i_con) = subplot(1,size(pl.mdata,2),i_contr);
    pl.clim = [-1 1] *max(abs(pl.mdata),[],'all');
    topoplot(pl.mdata(:,i_contr), EEG_cue.chanlocs, ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','off','colormap',flipud(cbrewer2('RdBu')),...
        'whitebk','on');
    title(sprintf('CDA | %s\n[%1.0f %1.0f]ms', ...
        pl.contrast_label{i_contr}, pl.time2plot), ...
         "Interpreter","none")
    colorbar
end

%% do some plotting CDA EEG all cons at one electrode cluster
pl.time2plot = [-100 1000];
pl.tidx = dsearchn(ERP_times', pl.time2plot');

pl.sub2plot = 1:numel(p.files2use);

% pl.elec2plot = {'CP3';'CP5';'P3';'P5';'P7';'P9';'PO3';'PO7'};
pl.elec2plot = {'CP4';'CP6';'P4';'P6';'P8';'P10';'PO4';'PO8'};
pl.elec2plot = {'PO4';'O2';'PO8';'I2'}; 
% pl.elec2plot = {'PO3';'O1';'PO7';'I1'};

pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({EEG_cue.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.data = squeeze(mean(ERP_data(pl.elec2plot_i,:,:,pl.sub2plot),[1,4]));

figure;
set(gcf,'Position',[100 100 900 500],'PaperPositionMode','auto')
plot(ERP_times,pl.data)
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
title("filtered ERP")
legend(p.con1label,"Interpreter","none","Location","Eastoutside")

h.a1 = axes('position',[0.74 0.68 0.25 0.25],'Visible','off');
topoplot(find(pl.elec2plot_i),EEG_cue.chanlocs,'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',4,1});

%% do some plotting CDA EEG contrast differences at one electrode cluster
pl.time2plot = [-100 1000];
pl.tidx = dsearchn(ERP_times', pl.time2plot');

pl.sub2plot = 1:numel(p.files2use);

pl.contrast = p.con3idx;

% pl.elec2plot = {'CP3';'CP5';'P3';'P5';'P7';'P9';'PO3';'PO7'};
pl.elec2plot = {'CP4';'CP6';'P4';'P6';'P8';'P10';'PO4';'PO8'};
pl.elec2plot = {'PO4';'O2';'PO8';'I2'}; 
pl.elec2plot = {'PO3';'O1';'PO7';'I1'};

pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({EEG_cue.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.data = nan([size(ERP_data, 2) size(pl.contrast,1) numel(pl.sub2plot)]);

for i_contr = 1:size(pl.contrast,1)
    pl.data(:,i_contr,:)=mean(ERP_data(pl.elec2plot_i,:,pl.contrast{i_contr}(1),:),1) - mean(ERP_data(pl.elec2plot_i,:,pl.contrast{i_contr}(2),:),1);
end

figure;
set(gcf,'Position',[100 100 900 500],'PaperPositionMode','auto')
plot(ERP_times,mean(pl.data,3))
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
title("filtered ERP | difference between cue left and right")
legend(p.con3label,"Interpreter","none","Location","Eastoutside")

h.a1 = axes('position',[0.74 0.68 0.25 0.25],'Visible','off');
topoplot(find(pl.elec2plot_i),EEG_cue.chanlocs,'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',4,1});

%% do some plotting CDA EEG contrast differences (cue left - cue right) electrode cluster (right - left)
pl.time2plot = [-100 1000];
pl.tidx = dsearchn(ERP_times', pl.time2plot');

pl.sub2plot = 1:numel(p.files2use);

pl.contrast = p.con3idx;

pl.elec2plot = {{'CP3';'CP5';'P3';'P5';'P7';'P9';'PO3';'PO7'}; {'CP4';'CP6';'P4';'P6';'P8';'P10';'PO4';'PO8'}};
pl.elec2plot = {{'PO3';'O1';'PO7';'I1'};{'PO4';'O2';'PO8';'I2'}};

pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({EEG_cue.chanlocs.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);

pl.data = nan([size(ERP_data, 2) size(pl.contrast,1) numel(pl.sub2plot)]);

for i_contr = 1:size(pl.contrast,1)
    pl.data(:,i_contr,:)= (mean(ERP_data(pl.elec2plot_i{2},:,pl.contrast{i_contr}(1),:),1) - mean(ERP_data(pl.elec2plot_i{2},:,pl.contrast{i_contr}(2),:),1)) - ...
        (mean(ERP_data(pl.elec2plot_i{1},:,pl.contrast{i_contr}(1),:),1) - mean(ERP_data(pl.elec2plot_i{1},:,pl.contrast{i_contr}(2),:),1));
end

figure;
set(gcf,'Position',[100 100 900 500],'PaperPositionMode','auto')
plot(ERP_times,mean(pl.data,3))
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
title("filtered ERP | difference cue left - right for right - left cluster")
legend(p.con3label,"Interpreter","none","Location","Eastoutside")
grid on

h.a1 = axes('position',[0.74 0.68 0.25 0.25],'Visible','off');
topoplot(find(pl.elec2plot_i{1}),EEG_cue.chanlocs,'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i{1}|pl.elec2plot_i{2}),'o','r',4,1});


%% do some plotting CDA EEG all cons | left cluster - right cluster
pl.time2plot = [-100 1000];
pl.tidx = dsearchn(ERP_times', pl.time2plot');

pl.sub2plot = 1:numel(p.files2use);

pl.elec2plot = {{'P9';'P7';'PO7'}; {'P8';'P10';'PO8'}};
pl.elec2plot = {{'CP3';'CP5';'P3';'P5';'P7';'P9';'PO3';'PO7'}; {'CP4';'CP6';'P4';'P6';'P8';'P10';'PO4';'PO8'}};
pl.elec2plot = {{'PO3';'O1';'PO7';'I1'}; {'PO4';'O2';'PO8';'I2'}};


pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({EEG_cue.chanlocs.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);

pl.data = squeeze(mean(mean(ERP_data(pl.elec2plot_i{1},:,:,pl.sub2plot),1)-mean(ERP_data(pl.elec2plot_i{2},:,:,pl.sub2plot),1),4));

figure;
set(gcf,'Position',[100 100 900 500],'PaperPositionMode','auto')
plot(ERP_times,pl.data)
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
title("filtered ERP | left cluster - right cluster")
legend(p.con1label,"Interpreter","none","Location","Eastoutside")
h.a1 = axes('position',[0.74 0.68 0.25 0.25],'Visible','off');
topoplot(find(pl.elec2plot_i{1}),EEG_cue.chanlocs,'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i{1}|pl.elec2plot_i{2}),'o','r',4,1});


%% do some plotting CDA sub contrasts | conta - ipsi
pl.time2plot = [-100 1000];
pl.tidx = dsearchn(ERP_times', pl.time2plot');

pl.sub2plot = 1:numel(p.files2use);

pl.elec2plot = {{'P9';'P7';'PO7'}; {'P8';'P10';'PO8'}};
% pl.elec2plot = {{'P7'}; {'P8'}};
% pl.elec2plot = {{'PO3';'O1';'I1'}; {'PO4';'O2';'I2'}};
pl.elec2plot = {{'CP3';'CP5';'P3';'P5';'P7';'P9';'PO3';'PO7';'O1';'I1'}; {'CP4';'CP6';'P4';'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'}};
% pl.elec2plot = {{'PO3';'O1';'PO7';'I1'};{'PO4';'O2';'PO8';'I2'}};



pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({EEG_cue.chanlocs.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);

pl.ContraIpsi_idx = [2 1; 1 2; 2 1; 1 2; 2 1; 1 2];
% pl.ContraIpsi_idx = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1];

pl.data = nan([size(ERP_data, [2 3]) numel(pl.sub2plot)]);

for i_con = 1:size(ERP_data,3)
    pl.data(:,i_con,:) =  squeeze(mean(ERP_data(pl.elec2plot_i{pl.ContraIpsi_idx(i_con,1)},:,i_con,pl.sub2plot),1) - ...
        mean(ERP_data(pl.elec2plot_i{pl.ContraIpsi_idx(i_con,2)},:,i_con,pl.sub2plot),1));
end

% average across sides
pl.data2 = nan(size(pl.data)./[1 2 1]);
for i_con = 1:size(p.con3idx,1)
    pl.data2(:,i_con,:) = mean(pl.data(:,p.con3idx{i_con},:),2);
end


figure;
set(gcf,'Position',[100 100 900 300],'PaperPositionMode','auto')
plot(ERP_times,mean(pl.data,3))
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
title("filtered ERP | contra - ipsi")
legend(p.con1label,"Interpreter","none","Location","Eastoutside")
grid on
h.a1 = axes('position',[0.74 0.68 0.25 0.25],'Visible','off');
topoplot(find(pl.elec2plot_i{1}),EEG_cue.chanlocs,'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i{1}|pl.elec2plot_i{2}),'o','r',4,1});


figure;
set(gcf,'Position',[100 100 500 200],'PaperPositionMode','auto')
plot(ERP_times,mean(pl.data2,3))
hold on; plot(ERP_times,mean(pl.data2,[2 3]),'Color',[0.2 0.2 0.2])
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
title("filtered ERP | contra - ipsi averaged across sides")
legend([p.con3label; {'mean'}],"Interpreter","none","Location","SouthEastoutside")
grid on
h.a1 = axes('position',[0.65 0.53 0.41 0.41],'Visible','off');
topoplot(find(pl.elec2plot_i{1}),EEG_cue.chanlocs,'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i{1}|pl.elec2plot_i{2}),'o','r',4,1});



%% extract some data | ERP | contra and ipsi
pl.time2plot = [-100 1000];
pl.tidx = dsearchn(ERP_times', pl.time2plot');

pl.sub2plot = 1:numel(p.files2use);

pl.elec2plot = {{'P9';'P7';'PO7'}; {'P8';'P10';'PO8'}};
% pl.elec2plot = {{'P7'}; {'P8'}};
% pl.elec2plot = {{'PO3';'O1';'I1'}; {'PO4';'O2';'I2'}};
pl.elec2plot = {{'CP3';'CP5';'P3';'P5';'P7';'P9';'PO3';'PO7';'O1';'I1'}; {'CP4';'CP6';'P4';'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'}};
% pl.elec2plot = {{'PO3';'O1';'PO7';'I1'};{'PO4';'O2';'PO8';'I2'}};

pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({EEG_cue.chanlocs.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);

% index contralateral electrode cluster
pl.ContraIpsi_idx = [2 1; 1 2; 2 1; 1 2; 2 1; 1 2];
% pl.ContraIpsi_idx = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1];
pl.contraIpsi_label = {'contra';'ipsi'};

pl.time2plot = [{[-1000 0]} {'tracking'};{[500 1000]} {'retention'};{[1000 1500]} {'probe'}];

% data: cluster x condition x time x participants
pl.data = nan(size(pl.elec2plot,1),size(ERP_data,3), size(pl.time2plot,1), numel(pl.sub2plot));

plout = [];
for i_elec = 1:size(pl.elec2plot,1)
    % index electrode cluster
    plout.elecs{i_elec} = vararg2str(pl.elec2plot{i_elec});
    for i_con =1:size(ERP_data,3)
        % index contra or ipsi
        plout.elec_label{i_elec,i_con} = pl.contraIpsi_label{pl.ContraIpsi_idx(i_con,i_elec)};
        % index conditions
        plout.con_idx(i_elec,i_con) = i_con;
        plout.con_label1{i_elec,i_con} = p.con1label{i_con};
        plout.con_label2{i_elec,i_con} = p.con2label{cell2mat(cellfun(@(x) any(ismember(x,i_con)), p.con2idx, 'UniformOutput', false))};
        plout.con_label3{i_elec,i_con} = p.con3label{cell2mat(cellfun(@(x) any(ismember(x,i_con)), p.con3idx, 'UniformOutput', false))};
        for i_time = 1:size(pl.time2plot,1)
            % index time windows
            pl.timeidx = dsearchn(ERP_times',pl.time2plot{i_time,1}');
            % write out time information
            plout.time{i_elec,i_con,i_time} = sprintf('[%1.0f %1.0f]ms',(pl.time2plot{i_time,1}));
            plout.time_label(i_elec,i_con, i_time) = pl.time2plot(i_time,2);
            % extract data
            pl.data(i_elec,i_con,i_time,:) = squeeze(mean( ...
                ERP_data(pl.elec2plot_i{i_elec},pl.timeidx(1):pl.timeidx(2),i_con,pl.sub2plot), ...
                [1,2]));
            % index participant
            plout.participant(i_elec,i_con, i_time,:) = p.files(p.files2use(pl.sub2plot));
        end
    end
end

% expand data
plout.elecs = repmat(plout.elecs', [1 size(pl.data,2:4)]);
plout.elec_label = repmat(plout.elec_label, [1 1 size(pl.data,3:4)]);
plout.con_idx = repmat(plout.con_idx, [1 1 size(pl.data,3:4)]);
plout.con_label1 = repmat(plout.con_label1, [1 1 size(pl.data,3:4)]);
plout.con_label2 = repmat(plout.con_label2, [1 1 size(pl.data,3:4)]);
plout.con_label3 = repmat(plout.con_label3, [1 1 size(pl.data,3:4)]);
plout.time = repmat(plout.time, [1 1 1 size(pl.data,4)]);
plout.time_label = repmat(plout.time_label, [1 1 size(pl.data,4)]);

plout.all_table = table(plout.participant(:), plout.elecs(:), plout.elec_label(:), plout.con_idx(:), ...
    plout.con_label1(:), plout.con_label2(:), plout.con_label3(:), plout.time(:), plout.time_label(:), pl.data(:), ...
    'VariableNames',{'participants','elecs','elec_label','con_idx','con_label1','con_label2','con_label3', ...
    'time','time_label','amplitude'});


t.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_workmem\data_in';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% write to textfile
writetable(plout.all_table,fullfile(t.path,sprintf('ERP_pilotdata_LargeLateralClust_%s.csv',t.datestr)),'Delimiter',';')


%% extract some data | SSVEP 
pl.sub2plot = 1:numel(p.files2use);

pl.elec2plot = {{'PO4';'PO8';'O2'}; {'P7';'PO7';'P9';'O1';'I1';'Iz'}};    % left stims one cluster and right stim 2 cluster

pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({EEG_cue.chanlocs.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);

pl.poslabel = {'left';'right'};
pl.con_posatt = [1 2 1 2 1 2];

pl.extRDKid= {[1 3 4 6]; [1 3 4 6]; [1 2 4 5]; [1 2 4 5]; [1 4]; [1 4]};
pl.extRDKposid = {[1 1 2 2];[1 1 2 2];[1 1 2 2];[1 1 2 2];[1 2];[1 2]};
pl.extRDKcols = {[1 2 1 2]; [1 2 1 2];[1 1 1 1];[1 1 1 1];[1 1];[1 1]};

pl.time2plot = [1 2];
pl.time2plot_label = {'preCue';'postCue'};

pl.freqrange = [-0.1 0.1];

pl.datastruct = [];
i_count = 1;

for i_sub = 1:numel(pl.sub2plot)
    for i_time = 1:numel(pl.time2plot)
        for i_con = 1:numel(p.con1idx)
            for i_rdk = 1:numel(pl.extRDKid{i_con})
                % write out participant
                pl.datastruct(i_count).participant = p.files{p.files2use(pl.sub2plot(i_sub))};
                % write out time
                pl.datastruct(i_count).time = sprintf('[%1.0f %1.0f]ms',p.cue_win2an{pl.time2plot(i_time)}*1000);
                pl.datastruct(i_count).time_label = pl.time2plot_label{pl.time2plot(i_time)};
                % write out condition
                pl.datastruct(i_count).con = p.con1idx(i_con);
                pl.datastruct(i_count).con_label1 = p.con1label{i_con};
                pl.datastruct(i_count).con_label2 = p.con2label{cell2mat(cellfun(@(x) any(ismember(x,i_con)), p.con2idx, 'UniformOutput', false))};
                pl.datastruct(i_count).con_label3 = p.con3label{cell2mat(cellfun(@(x) any(ismember(x,i_con)), p.con3idx, 'UniformOutput', false))};
                % write RDK infos
                pl.datastruct(i_count).RDK_id = pl.extRDKid{i_con}(i_rdk);
                pl.datastruct(i_count).RDK_posid = pl.extRDKposid{i_con}(i_rdk);
                pl.datastruct(i_count).RDK_poslabel = pl.poslabel{pl.datastruct(i_count).RDK_posid};
                pl.datastruct(i_count).RDK_freq = behavior.RDK.RDK(pl.datastruct(i_count).RDK_id).freq;
                % position attended
                if pl.extRDKposid{pl.datastruct(i_count).con}(i_rdk)==pl.con_posatt(pl.datastruct(i_count).con)
                    pl.datastruct(i_count).RDK_posatt = 'attended';
                else
                    pl.datastruct(i_count).RDK_posatt = 'unattended';
                end
                % color attended
                if pl.extRDKcols{pl.datastruct(i_count).con}(i_rdk)==1
                    pl.datastruct(i_count).RDK_colatt = 'attended';
                else
                    pl.datastruct(i_count).RDK_colatt = 'unattended';
                end
                % extract electrodes
                pl.datastruct(i_count).RDK_elec = vararg2str(pl.elec2plot{pl.datastruct(i_count).RDK_posid});
                % extract amplitude
                freqidx = dsearchn(FFT_freqs', (pl.freqrange+pl.datastruct(i_count).RDK_freq)');
                pl.datastruct(i_count).RDK_amp = squeeze(mean(FFT_data_evo( ...
                    freqidx(1):freqidx(2), ...
                    pl.elec2plot_i{pl.datastruct(i_count).RDK_posid}, ...
                    pl.datastruct(i_count).con, ...
                    pl.time2plot(i_time), ...
                    pl.sub2plot(i_sub)), ...
                    [1 2]));
                i_count = i_count+1;

            end
        end
    end
end

pl.datatable = struct2table(pl.datastruct);


t.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_workmem\data_in';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% write to textfile
writetable(pl.datatable,fullfile(t.path,sprintf('FFT_SSVEP_pilotdata_%s.csv',t.datestr)),'Delimiter',';')
