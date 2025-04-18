%% script to analyze behavioral pilot data for ssvep_workmem
clearvars

p.path =                'N:\AllgPsy\experimental_data\2025_SSVEP_workmem\behavior\';
p.path =                '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2025_SSVEP_workmem\behavior\';

p.files =               {'Pilot_CG_VP90_timing.mat';'Pilot_Chris_VP92_timing.mat';'Pilot_Sabrina_VP96_timing.mat';
                            'Pilot_Katha_VP94_timing.mat';'Pilot_Christian_VP97_timing.mat';'VP01_timing.mat';'VP02_timing.mat';'VP03_timing.mat'};

p.responsewin_main =    [0.2 1.5];
p.responsewin_pre =     [0.2 1]; % according to p.targ_respwin from run_posnalpha


p.con1idx =             [1 2 3 4 5 6];    
p.con1label =           {'cue_left_2t2d'; 'cue_right_2t2d'; 'cue_left_4t'; 'cue_right_4t'; 'cue_left_2t'; 'cue_right_2t'};
                        % [RDK1 RDK3; RDK4 RDK6] [RDK1 RDK3; RDK4 RDK6] [RDK1 RDK2; RDK4 RDK5] [RDK1 RDK2; RDK4 RDK5] [RDK1; RDK4] [RDK1; RDK4] 
p.con2idx =             {[1 3 5];[2 4 6]};
p.con2label =           {'cue_left';'cue_right'};
p.con3idx =             {[1 2];[3 4]; [5 6]};
p.con3label =           {'2t_2d';'4t';'2t'};

p.colnum =              [0.4706; 0]; % num of red channel
p.collabel =            {'orange';'blue'};


%% actual calculation
for i_sub = 1:numel(p.files)
    i_sub
    % load data
    data_in.resp.experiment = repmat({[nan]},1,8);
    data_in.button_presses.experiment = repmat({[nan]},1,8);

    temp.data_in=load(sprintf('%s%s',p.path,p.files{i_sub}));
    data_in.RDK =  temp.data_in.RDK.RDK;
    
    temp.index=cell2mat(cellfun(@(x) ~isempty(cell2mat({x(:).trialnumber})), temp.data_in.resp.experiment,'UniformOutput',false));
    data_in.resp.experiment(temp.index)=temp.data_in.resp.experiment(temp.index);
    data_in.button_presses.experiment(temp.index)=temp.data_in.button_presses.experiment(temp.index);
    
       
    %% setup response matrices
    data_all_struct = horzcat(data_in.resp.experiment{:});
    data_all_table = struct2table(data_all_struct);
    clear data_events
    
    % add data
    data_all_table.participant = repmat({p.files{i_sub}(1:end-11)},size(data_all_table,1),1);

    temp.col = p.collabel(ismember(p.colnum,data_in.RDK(1).col(1)));
    data_all_table.relevantcol = repmat(temp.col,size(data_all_table,1),1);
    
    % write into one file   
    if i_sub == 1
        response = data_all_table;
    else
        response = vertcat(response,data_all_table);
    end
end

%% check for some measures
% interim statistics
clear interim
for i_sub = 1:numel(p.files)
    % index participant
    temp.idx=strcmp(response.participant,p.files{i_sub}(1:end-11));

    % index participant
    interim(i_sub).participant = {p.files{i_sub}(1:end-11)};

    % extract hitrate
    interim(i_sub).hitrate_all = sum(strcmp(response.event_response_type(temp.idx),'hit'))/sum(temp.idx);
    interim(i_sub).errorrate_all = sum(strcmp(response.event_response_type(temp.idx),'error'))/sum(temp.idx);
    % for conditions
    for i_con = 1:numel(p.con3label)
        % index condition
        temp.idx2 = ismember(response.condition,p.con3idx{i_con});
%         temp.name = sprintf('hitrate_%s',p.con3label{i_con});
%         interim(i_sub) = setfield(interim(i_sub), temp.name, sum(strcmp(response.event_response_type(temp.idx & temp.idx2), 'hit')) / sum(temp.idx & temp.idx2));
        interim(i_sub).(sprintf('hitrate_%s',p.con3label{i_con})) = ...
            sum(strcmp(response.event_response_type(temp.idx&temp.idx2),'hit'))/sum(temp.idx&temp.idx2);
        interim(i_sub).(sprintf('errorrate_%s',p.con3label{i_con})) = ...
            sum(strcmp(response.event_response_type(temp.idx&temp.idx2),'error'))/sum(temp.idx&temp.idx2);
    end
    
    % RT
    interim(i_sub).RT_hit_mean = mean(response.event_response_RT(strcmp(response.event_response_type,'hit')&temp.idx));
    interim(i_sub).RT_hit_std = std(response.event_response_RT(strcmp(response.event_response_type,'hit')&temp.idx));
    interim(i_sub).RT_error_mean = mean(response.event_response_RT(strcmp(response.event_response_type,'error')&temp.idx));
    interim(i_sub).RT_error_std = std(response.event_response_RT(strcmp(response.event_response_type,'error')&temp.idx));
    for i_con = 1:numel(p.con3label)
        % index condition
        temp.idx2 = ismember(response.condition,p.con3idx{i_con});
        interim(i_sub).(sprintf('RT_hit_mean_%s',p.con3label{i_con})) = ...
            mean(response.event_response_RT(strcmp(response.event_response_type,'hit')&temp.idx&temp.idx2));
        interim(i_sub).(sprintf('RT_hit_std_%s',p.con3label{i_con})) = ...
            std(response.event_response_RT(strcmp(response.event_response_type,'hit')&temp.idx&temp.idx2));
        interim(i_sub).(sprintf('RT_error_mean_%s',p.con3label{i_con})) = ...
            mean(response.event_response_RT(strcmp(response.event_response_type,'error')&temp.idx&temp.idx2));
        interim(i_sub).(sprintf('RT_error_std_%s',p.con3label{i_con})) = ...
            std(response.event_response_RT(strcmp(response.event_response_type,'error')&temp.idx&temp.idx2));
    end
end




%% summary statistics/descriptives
% remove data that is not relevant
response_out = removevars(response, {'color_attended';'eventcolor';'stimangles';'eventangles';'button_presses_fr';'button_presses_t'});
% adapt certain columm
response_out.RDK2display = cellfun(@(x) vararg2str(x),response_out.RDK2display,'UniformOutput',false);
response_out.RDK2display = cellfun(@(x) x(1:end-2),response_out.RDK2display,'UniformOutput',false);
response_out.con1label = arrayfun(@(x) p.con1label(x),response_out.condition);
response_out.con2label = arrayfun(@(x) p.con2label(any(ismember(cell2mat(p.con2idx),x),2)), response_out.condition);
response_out.con3label = arrayfun(@(x) p.con3label(any(ismember(cell2mat(p.con3idx),x),2)), response_out.condition);



% general summary mean RT and hit rate
% export data
p.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_workmem\data_in';
p.filename = 'behavior.csv';
writetable(response_out, fullfile(p.path,p.filename))

