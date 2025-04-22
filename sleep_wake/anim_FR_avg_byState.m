% go through multi animal cell structure with state data and create animal
% average FR plots colored by behavioral state

% add extwake to sleep transition avg FR and extsleep to wake 

cd('Z:\BrianCary\CODE\Dan_SW_FR_analysis\Dan_version_ofAlejCode')

filepath = "\\files.brandeis.edu\turrigiano-lab\DPL\PROJECTS\StateCoding\MasterStrct_StableBase_Hr66_withStates_v8_BL71incl.mat";
filepath = "Z:\DPL\PROJECTS\StateCoding\MasterStrct_StableBase_Hr66_withStates_v8_BL71incl.mat";

% filepath = "Z:\BrianCary\CODE\Dan_SW_FR_analysis\Dan_version_ofAlejCode\stored_data\ContCell_MasterStrct_v3_DL78pred.mat";
% 
% filepath = "Z:\DPL\PROJECTS\StateCoding\STATEvsMOD\PosModStableBase_Hr66.mat";
% filepath = "Z:\DPL\PROJECTS\StateCoding\STATEvsMOD\NegModStableBase_Hr66.mat";
% filepath = "Z:\DPL\PROJECTS\StateCoding\STATEvsMOD\pFS\pFSStableBase_Hr66.mat";


loadFile = load(filepath);

loadFields = fieldnames(loadFile);

MASTER = loadFile.(loadFields{1}).MASTER;
try
    STATES = loadFile.(loadFields{1}).STATETIMES;
catch
    STATES = loadFile.(loadFields{1}).STATES;
end
DAYSTART = STATES.DAYSTART;
anim_fields = fieldnames(STATES);
anim_fields(strcmp(anim_fields,'DAYSTART')) = [];
n_anims = numel(anim_fields); % minus daystart field

cell_anims = {MASTER.animal};

%% initialize vars
PLOT_EACH_CELL = 1;
PLOT_EACH_ANIM = 0;

setFigureDefaults

% colors for plotting
%% 
c_rem   = [25 181 149]./255;
c_nrem  = [131 49 146]./255;
c_aw    = [201 28 101]./255;
c_qw    = [247 148 41]./255;

c_rem_cont   = c_rem*0.75;
c_nrem_cont  = c_nrem*0.75;
c_aw_cont    = c_aw*0.75;
c_qw_cont    = c_qw*0.75;

lightOn_hr = 7.5;

test_type = 'sign';

% Vid params
fr_int = 0.25; % in seconds

smooth_ON = 1;
gauss_wind_sec = 15;
win_size = round(gauss_wind_sec/fr_int); % make x sec gauss wind
gauss_wind = gausswin(win_size);
gauss_wind = gauss_wind/sum(gauss_wind);

baseline_hrs = [60, 72]; % 
anim_smooth_sec = 60;

% posthunt_mean_hrs = [72+8+28+12, 72+36+32];
posthunt_mean_hrs = [72+8+4+48, 72+36+36];
% posthunt_mean_hrs = [72+8+4, 72+36+36];

ext_wake_time_thresh = 30; % minutes
ext_sleep_time_thresh = 30; % minutes
max_state_interrupt = 120; %seconds


% initialize variables
cell_ct = 0;
REM_meanFRbyepoch = []; NREM_meanFRbyepoch = []; AW_meanFRbyepoch = []; QW_meanFRbyepoch = [];
REM_meanCVbyepoch = []; NREM_meanCVbyepoch = []; AW_meanCVbyepoch = []; QW_meanCVbyepoch = [];
REM_meanDUR_byepoch = []; NREM_meanDUR_byepoch = []; AW_meanDUR_byepoch = []; QW_meanDUR_byepoch = [];
REM_meanFRbycell_BL = [];
NREM_meanFRbycell_BL = [];
AW_meanFRbycell_BL = [];
QW_meanFRbycell_BL = [];

REM_meanCVbycell_BL = [];
NREM_meanCVbycell_BL = [];
AW_meanCVbycell_BL = [];
QW_meanCVbycell_BL = [];

REM_meanFRbycell_hunt = [];
NREM_meanFRbycell_hunt = [];
AW_meanFRbycell_hunt = [];
QW_meanFRbycell_hunt = [];

REM_meanCVbycell_hunt = [];
NREM_meanCVbycell_hunt = [];
AW_meanCVbycell_hunt = [];
QW_meanCVbycell_hunt = [];    

for anim_i = 1:n_anims
    %% init
    anim_states = STATES.(anim_fields{anim_i});
    anim_daystart = DAYSTART.(anim_fields{anim_i});
%     state_start = anim_states(1,2);

    anim_inds = strcmp(cell_anims,anim_fields{anim_i});
    anim_cells = MASTER(anim_inds);

%     if isempty(anim_cells)
%         continue
%     end

    expt_start = anim_cells(1).EXPTSTART;
    trem = anim_cells(1).trem;
      
    last_state_time = anim_states(end,2) - expt_start;
    time_bins = (trem:fr_int:last_state_time) + trem + anim_daystart*24*3600;

    % create state vec
    %1 = REM, 2 = NREM, 4 = AW, 5 = QW
    state_times = anim_states(:,2);
    state_times = state_times - expt_start + trem + anim_daystart*24*3600;
    state_vec = nan(1,length(time_bins));
    for s = 1:length(state_times)-1
        s_inds = time_bins > state_times(s) & time_bins <= state_times(s+1);
        state_vec(s_inds) = anim_states(s,1);
    end
    state_vec = state_vec(2:end); % remove first element to account for binning
    

    %% cell loop
    n_cells = length(anim_cells);
    cell_fr_mat = nan([n_cells, length(time_bins)-1]);
    cell_fr_mat_sm = nan([n_cells, length(time_bins)-1]);
    for cell_i = 1:n_cells
        
        this_cell = anim_cells(cell_i);
        
        all_sp_times = this_cell.time + anim_daystart*24*3600;
        onTimes = this_cell.onTime*3600 + anim_daystart*24*3600;
        offTimes = this_cell.offTime*3600 + anim_daystart*24*3600;
        
        % re-center time bins
        time_bins_sh = time_bins(2:end) - fr_int/2;

        cell_sp_times = [];
        num_ontimes = length(onTimes);
        nonnan_inds = [];
        for ii = 1:num_ontimes
            add_sps = all_sp_times(all_sp_times > onTimes(ii)...
                & all_sp_times < offTimes(ii));
            cell_sp_times = [cell_sp_times; add_sps];
            
            nonnan_ind_ii = find(time_bins_sh > onTimes(ii) & time_bins_sh <  offTimes(ii))';
            nonnan_inds = [nonnan_inds; nonnan_ind_ii];
        end       

%         if isempty(time_bins)
%             continue
%         end
%         cell_sp_times = cell_sp_times + anim_daystart*24*3600;
        cell_fr = histcounts(cell_sp_times,time_bins) ./ fr_int;

        nonnan_vec = logical(zeros([1,length(cell_fr)]));
        nonnan_vec(nonnan_inds) = 1;
        cell_fr(~nonnan_vec) = NaN;

        % apply gauss filt to cell FR
        if smooth_ON == 1
            cell_fr = filter(gauss_wind,1,cell_fr);
        end 
        
        this_cell.time = this_cell.time + anim_daystart*24*3600;
        this_cell.onTime = this_cell.onTime*3600 + anim_daystart*24*3600;
        this_cell.offTime = this_cell.offTime*3600 + anim_daystart*24*3600;
        cell_statetimes = adjust_statetimes_recov(anim_states,this_cell,anim_daystart);
        use_daystart = 0;
        nsec = 0;
        [REM_mean_cell, NREM_mean_cell, AW_mean_cell, QW_mean_cell,...
            REM_mean_epochs, NREM_mean_epochs, AW_mean_epochs, QW_mean_epochs, ...
            REM_cvcell, NREM_cvcell, AW_cvcell, QW_cvcell,...
            REM_cvepoch, NREM_cvepoch, AW_cvepoch, QW_cvepoch,...
            REM_durcell, NREM_durcell, AW_durcell, QW_durcell,...
            REM_durepoch, NREM_durepoch, AW_durepoch, QW_durepoch, zerocount] = ...
            SW_meanFR_4state(this_cell, cell_statetimes, nsec, use_daystart);
        
        % baseline
        this_cell_BL = this_cell;
        non_bl_inds = (this_cell_BL.time < baseline_hrs(1)*3600) | (this_cell_BL.time > baseline_hrs(2)*3600);
        this_cell_BL.time(non_bl_inds) = [];
        use_daystart = 0;
        nsec = 0;
        [REM_mean_cell_BL, NREM_mean_cell_BL, AW_mean_cell_BL, QW_mean_cell_BL,...
            ~, ~, ~, ~, ...
            REM_cvcell_BL, NREM_cvcell_BL, AW_cvcell_BL, QW_cvcell_BL,...
            ~, ~, ~, ~,...
            ~, ~, ~, ~,...
            ~, ~, ~, ~, ~] = ...
            SW_meanFR_4state(this_cell_BL, cell_statetimes, nsec, use_daystart);

        % hunt
        this_cell_hunt = this_cell;
        non_hunt_inds = (this_cell_hunt.time < posthunt_mean_hrs(1)*3600) | (this_cell_hunt.time > posthunt_mean_hrs(2)*3600);
        this_cell_hunt.time(non_hunt_inds) = [];
        use_daystart = 0;
        nsec = 0;
        [REM_mean_cell_hunt, NREM_mean_cell_hunt, AW_mean_cell_hunt, QW_mean_cell_hunt,...
            ~, ~, ~, ~, ...
            REM_cvcell_hunt, NREM_cvcell_hunt, AW_cvcell_hunt, QW_cvcell_hunt,...
            ~, ~, ~, ~,...
            ~, ~, ~, ~,...
            ~, ~, ~, ~, ~] = ...
            SW_meanFR_4state(this_cell_hunt, cell_statetimes, nsec, use_daystart);

        % FR
        cell_ct = cell_ct + 1;
%         cell_ct = cell_i; % if you want to refactor for separate animal
%         plotting

        REM_meanFRbycell(cell_ct,1)   = REM_mean_cell;
        REM_meanFRbyepoch        = [REM_meanFRbyepoch; REM_mean_epochs];
        NREM_meanFRbycell(cell_ct,1)  = NREM_mean_cell;
        NREM_meanFRbyepoch       = [NREM_meanFRbyepoch; NREM_mean_epochs];
        AW_meanFRbycell(cell_ct,1)    = AW_mean_cell;
        AW_meanFRbyepoch         = [AW_meanFRbyepoch; AW_mean_epochs];
        QW_meanFRbycell(cell_ct,1)    = QW_mean_cell;
        QW_meanFRbyepoch         = [QW_meanFRbyepoch; QW_mean_epochs];
        
        
        % CV
        REM_meanCVbycell(cell_ct,1)   = REM_cvcell;
        REM_meanCVbyepoch        = [REM_meanCVbyepoch; REM_cvepoch];
        NREM_meanCVbycell(cell_ct,1)  = NREM_cvcell;
        NREM_meanCVbyepoch       = [NREM_meanCVbyepoch; NREM_cvepoch];
        AW_meanCVbycell(cell_ct,1)    = AW_cvcell;
        AW_meanCVbyepoch         = [AW_meanCVbyepoch; AW_cvepoch];
        QW_meanCVbycell(cell_ct,1)    = QW_cvcell;
        QW_meanCVbyepoch         = [QW_meanCVbyepoch; QW_cvepoch];
        
        % EPOCH DURATION
        REM_meanDUR_bycell(cell_ct,1)   = REM_durcell;
        REM_meanDUR_byepoch        = [REM_meanDUR_byepoch; REM_durepoch];
        NREM_meanDUR_bycell(cell_ct,1)  = NREM_durcell;
        NREM_meanDUR_byepoch       = [NREM_meanDUR_byepoch; NREM_durepoch];
        AW_meanDUR_bycell(cell_ct,1)    = AW_durcell;
        AW_meanDUR_byepoch         = [AW_meanDUR_byepoch; AW_durepoch];
        QW_meanDUR_bycell(cell_ct,1)    = QW_durcell;
        QW_meanDUR_byepoch         = [QW_meanDUR_byepoch; QW_durepoch];


        REM_meanFRbycell_BL(cell_ct,1)   = REM_mean_cell_BL;
        NREM_meanFRbycell_BL(cell_ct,1)  = NREM_mean_cell_BL;
        AW_meanFRbycell_BL(cell_ct,1)    = AW_mean_cell_BL;
        QW_meanFRbycell_BL(cell_ct,1)    = QW_mean_cell_BL;
        
        REM_meanCVbycell_BL(cell_ct,1)   = REM_cvcell_BL;
        NREM_meanCVbycell_BL(cell_ct,1)  = NREM_cvcell_BL;
        AW_meanCVbycell_BL(cell_ct,1)    = AW_cvcell_BL;
        QW_meanCVbycell_BL(cell_ct,1)    = QW_cvcell_BL;

        REM_meanFRbycell_hunt(cell_ct,1)   = REM_mean_cell_hunt;
        NREM_meanFRbycell_hunt(cell_ct,1)  = NREM_mean_cell_hunt;
        AW_meanFRbycell_hunt(cell_ct,1)    = AW_mean_cell_hunt;
        QW_meanFRbycell_hunt(cell_ct,1)    = QW_mean_cell_hunt;
        
        REM_meanCVbycell_hunt(cell_ct,1)   = REM_cvcell_hunt;
        NREM_meanCVbycell_hunt(cell_ct,1)  = NREM_cvcell_hunt;
        AW_meanCVbycell_hunt(cell_ct,1)    = AW_cvcell_hunt;
        QW_meanCVbycell_hunt(cell_ct,1)    = QW_cvcell_hunt;        

        if PLOT_EACH_CELL == 1
            % make figure
            fig_lines1 = figure('Position',[0.0271 0.4704 0.5469 0.4315]);
            ax_lines1 = axes('Parent',fig_lines1);
            hold(ax_lines1,'on')
    
            % plot initial lines
            add_LD_bars(ax_lines1,12,6,[0 time_bins_sh(end)/3600+12]);
        
            plot(ax_lines1,time_bins_sh./3600,cell_fr,'color',[0 0 0 0.7]);       
            ylabel(ax_lines1,'Firing Rate (Hz)')
            xlabel(ax_lines1,'Time (hrs)')

            nrem_fr = cell_fr;
            rem_fr = cell_fr;
            qw_fr = cell_fr;
            aw_fr = cell_fr;
            nrem_fr(state_vec ~= 2) = NaN;
            rem_fr(state_vec ~= 1) = NaN;
            qw_fr(state_vec ~= 5) = NaN;
            aw_fr(state_vec ~= 4) = NaN;
            
            line_objs1= plot(ax_lines1,time_bins_sh./3600,nrem_fr,'color',c_nrem);
            line_objs2 = plot(ax_lines1,time_bins_sh./3600,rem_fr,'color',c_rem);  
            line_objs3 = plot(ax_lines1,time_bins_sh./3600,qw_fr,'color',c_qw);
            line_objs4 = plot(ax_lines1,time_bins_sh./3600,aw_fr,'color',c_aw);

            title_str = ['Cell: ',num2str(cell_i)];
            title(ax_lines1,title_str)        
                
            box(ax_lines1,'off')
            leg_names = {'nrem','rem','qw','aw'};
            legend(ax_lines1,[line_objs1,line_objs2,line_objs3,line_objs4],leg_names);
            
            xlim(ax_lines1, [0, time_bins_sh(end)/3600+12])

%             xlim([2.2454e5    2.3378e5]./3600)

            ylim(ax_lines1, [0, max(cell_fr)*1.2])

        end

        baseline_inds =  time_bins_sh > baseline_hrs(1)*3600 & time_bins_sh <=  baseline_hrs(2)*3600;
        norm_cell_fr = cell_fr/mean(cell_fr(baseline_inds),'omitnan');
        cell_fr_mat(cell_i,:) = norm_cell_fr;

        nan_inds = isnan(norm_cell_fr);
        norm_cell_fr_nonan = norm_cell_fr;
        norm_cell_fr_nonan(nan_inds) = 0;
        norm_cell_fr_sm = smooth(norm_cell_fr_nonan, anim_smooth_sec*(1/fr_int));
        norm_cell_fr_sm(nan_inds) = NaN;
        cell_fr_mat_sm(cell_i,:) = norm_cell_fr_sm;
    end
    
    %% plotting
    if PLOT_EACH_ANIM == 1
        anim_fr = nanmean(cell_fr_mat_sm,1);
    
    %     if isempty(anim_fr)
    %         continue
    %     end
    
        cell_fr_mat_sm(cell_fr_mat_sm==0) = NaN;
        
    %     % optional removal of specific sections
    %     t0 = 85.022;
    %     t1 = 85.084;
    %     inds_to_remove = time_bins_sh/3600 > t0 & time_bins_sh/3600 < t1;
    %     cell_fr_mat_sm(:,inds_to_remove) = NaN;
    
    
        anim_nrem_fr = cell_fr_mat_sm;
        anim_rem_fr = cell_fr_mat_sm;
        anim_qw_fr = cell_fr_mat_sm;
        anim_aw_fr = cell_fr_mat_sm;
        anim_nrem_fr(:,state_vec ~= 2) = NaN;
        anim_rem_fr(:,state_vec ~= 1) = NaN;
        anim_qw_fr(:,state_vec ~= 5) = NaN;
        anim_aw_fr(:,state_vec ~= 4) = NaN;
        
        adj_anim_states = [anim_states(:,1), state_times];
    
        % remove repeats
        st_diff = diff(adj_anim_states(:,1));
        kill_these = find(st_diff==0);
        adj_anim_states(kill_these+1,:) = [];
        
        % remove short states
        interrupt_threshold = 10;
        st_timediff = diff(adj_anim_states(:,2));
        too_short = find(st_timediff <= interrupt_threshold);
        adj_anim_states(too_short,:) = [];
        
        % remove repeats again
        st_diff = diff(adj_anim_states(:,1));
        kill_these = find(st_diff==0);
        adj_anim_states(kill_these+1,:) = [];
    
        ext_wakes = find_extended_wake(adj_anim_states,ext_wake_time_thresh,max_state_interrupt);
        ext_sleeps = find_extended_sleep(adj_anim_states,ext_sleep_time_thresh,max_state_interrupt);
    
        wake_durs = (adj_anim_states(ext_wakes(:,2),2) - adj_anim_states(ext_wakes(:,1),2));
        sleep_durs = (adj_anim_states(ext_sleeps(:,2),2) - adj_anim_states(ext_sleeps(:,1),2));
    
        time_bins_toplot = time_bins_sh - (baseline_hrs(2)+8)*3600;
    %     time_bins_toplot = time_bins_sh;
    
        % make figure
        fig_lines1 = figure('Position',[0.0271 0.4704 0.5469 0.4315]);
    %     fig_lines1 = figure('Position',[154 515 1195 458]);
        ax_lines1 = axes('Parent',fig_lines1);
        set(ax_lines1,'FontSize',18)
        hold(ax_lines1,'on')
    
        % plot initial lines
        add_LD_bars(ax_lines1,12 - (baseline_hrs(2)+8),7,[0 time_bins_toplot(end)/3600+12]);
    
    %     plot(ax_lines1,time_bins_sh./3600,cell_fr,'color',[0 0 0 0.7]);       
        ylabel(ax_lines1,'Firing Rate (Norm.)','FontSize',18)
        xlabel(ax_lines1,'Time Post Hunt (hrs)','FontSize',18)
        
        alpha = 0.75;
        line_objs1= plot(ax_lines1,time_bins_toplot./3600,mean(anim_nrem_fr,1,'omitnan'),'color',[c_nrem,alpha],'LineWidth',1.5);
        line_objs2 = plot(ax_lines1,time_bins_toplot./3600,mean(anim_rem_fr,1,'omitnan'),'color',[c_rem,alpha],'LineWidth',1.5);
        line_objs3 = plot(ax_lines1,time_bins_toplot./3600,mean(anim_qw_fr,1,'omitnan'),'color',[c_qw,alpha],'LineWidth',1.5);
        line_objs4 = plot(ax_lines1,time_bins_toplot./3600,mean(anim_aw_fr,1,'omitnan'),'color',[c_aw,alpha],'LineWidth',1.5);
    
        title_str = ['Anim: ',anim_fields{anim_i}];
        title(ax_lines1,title_str)        
            
        box(ax_lines1,'off')
        
        ymax = max(mean(anim_aw_fr,1,'omitnan'))*1.2;
        plot([0, 0],[0, ymax],'r--','LineWidth',2)
    
        leg_names = {'nrem','rem','quiet wake','active wake'};
        legend(ax_lines1,[line_objs1,line_objs2,line_objs3,line_objs4],leg_names,'FontSize',14);    
        xlim(ax_lines1, [-(baseline_hrs(2)+8), time_bins_toplot(end)/3600+12])
        ylim(ax_lines1, [0, ymax])
    
    %     xlim([15, 28.5])
    %     ylim([-0.2, 4])
    
    %     xlim([39, 53])
    %     ylim([-0.2 7])
        % draw extended wake boxes
        for e_i = 1:size(ext_wakes,1)
            
            e_i_st = adj_anim_states(ext_wakes(e_i,1),2)/3600- (baseline_hrs(2)+8);
            e_i_end = adj_anim_states(ext_wakes(e_i,2)+1,2)/3600- (baseline_hrs(2)+8);
            width = e_i_end - e_i_st;
            height = ymax*0.75;
    %         height = 4.5;
            rect_obj = rectangle(ax_lines1, 'Position', [e_i_st 0 width height]);
    
            rect_obj.LineWidth = 1;
            rect_obj.LineStyle = '--';
            rect_obj.EdgeColor = c_qw;
    
        end
    
        % draw extended sleep boxes
        for e_i = 1:size(ext_sleeps,1)
            
            e_i_st = adj_anim_states(ext_sleeps(e_i,1),2)/3600- (baseline_hrs(2)+8);
            e_i_end = adj_anim_states(ext_sleeps(e_i,2)+1,2)/3600- (baseline_hrs(2)+8);
            width = e_i_end - e_i_st;
            height = ymax*0.75;
    %         height = 3.25;
    
            rect_obj = rectangle(ax_lines1, 'Position', [e_i_st 0 width height]);
    
            rect_obj.LineWidth = 1;
            rect_obj.LineStyle = '--';
            rect_obj.EdgeColor = c_nrem;
    
        end
        
        
        wind_sec = 10*60;
        norm_sec = 2*60;
        edges_toplot = (fr_int:fr_int:wind_sec) - wind_sec/2;
        sleep2extwake = [];
        extwake2sleep = [];
        for e_i = 1:size(ext_wakes,1)
            
            e_i_st = adj_anim_states(ext_wakes(e_i,1),2);
            e_i_end = adj_anim_states(ext_wakes(e_i,2)+1,2);
            
            if e_i_st/3600 > baseline_hrs(1)
                sleep2wake_inds = find(time_bins_sh >= (e_i_st - wind_sec/2) & time_bins_sh < (e_i_st + wind_sec/2));
                wake2sleep_inds = find(time_bins_sh >= (e_i_end - wind_sec/2) & time_bins_sh < (e_i_end + wind_sec/2));
                
                try
                    if length(sleep2wake_inds) == length(edges_toplot)
    
                        % skip if zeros in anim_fr
                        if sum(anim_fr(sleep2wake_inds)==0) == 0 &&...
                                sum(anim_fr(wake2sleep_inds)==0) == 0
    
                            sleep2extwake(e_i,:) = anim_fr(sleep2wake_inds)...
                                    / nanmean(anim_fr(sleep2wake_inds(1:round(norm_sec/fr_int))));
                            extwake2sleep(e_i,:) = anim_fr(wake2sleep_inds)...
                                    / nanmean(anim_fr(wake2sleep_inds(1:round(norm_sec/fr_int))));
                        end
                    end
                end
            end
        end
        sleep2extwake(sleep2extwake==0) = NaN;
        extwake2sleep(extwake2sleep==0) = NaN;
    
        figure; hold on
        alpha = 0.3;
        sh_color = 'g';
        s1 = stdshade(sleep2extwake,alpha,sh_color,edges_toplot,1,'sem');    
        plot([0,0],[0.5, 2],'k')
        ylim([0.75, 1.5])
        ylabel('Norm FR')
        xlabel('Seconds since extwake start')
        title(['sleep to extwake transition FR - norm to sleep - anim=',anim_fields{anim_i}])
    
        figure; hold on
        alpha = 0.3;
        sh_color = 'g';
        s1 = stdshade(extwake2sleep,alpha,sh_color,edges_toplot,1,'sem');    
        plot([0,0],[0.5, 2],'k')
        ylim([0.5, 1.25])
        ylabel('Norm FR')
        xlabel('Seconds since extwake end')
        title(['extwake to sleep transition FR - norm to extwake - anim=',anim_fields{anim_i}])
    
    % 
    %     % plot average FR and CV in the 4 states for each cell
    %     figure('Position',[0.0646 0.2611 0.2323 0.5981])
    %     
    %     data_toplot = [AW_meanFRbycell,...
    %         QW_meanFRbycell,...
    %         NREM_meanFRbycell,...
    %         REM_meanFRbycell];
    % 
    %     Colors = [c_aw;c_qw;c_nrem; c_rem];
    %     [s_obj, x, y] = UnivarScatter_ATP(data_toplot,'BoxType','SEM',...
    %         'Width',1,'Compression',500,'MarkerFaceColor',Colors,...
    %         'PointSize',55,'StdColor','none',...
    %         'Whiskers','lines','WhiskerLineWidth',2.5,'MarkerEdgeColor',Colors);
    %     box off
    %     hold on
    %     plot(x', y','Color',[0 0 0 0.1],'LineWidth',2)
    %     
    %     xlim([0.5 4.5])
    %     ylim([0.03 10])
    %     box off
    %     set(gca,'yscale','log')
    %     ylabel('Firing Rate (Hz)')
    %     set(gca,'YTickLabel',[0.1, 1, 10])
    %     set(gca,'XTick',1:4)
    %     set(gca,'XTickLabel',{'Active','Quiet','NREM','REM'},'XTickLabelRotation',45)
    %     title(['Cell Avg. Fr by State -- Anim:',anim_fields{anim_i}],'FontSize',16)    
    
        drawnow

    end
    %%
end

%% full exp avgs
% plot average FR in the 4 states for each cell
figure('Position',[0.0646 0.2611 0.2323 0.5981])

data_toplot = [AW_meanFRbycell,...
    QW_meanFRbycell,...
    NREM_meanFRbycell,...
    REM_meanFRbycell];

Colors = [c_aw;c_qw;c_nrem; c_rem];
[s_obj, x, y] = UnivarScatter_ATP(data_toplot,'BoxType','SEM',...
    'Width',1,'Compression',500,'MarkerFaceColor',Colors,...
    'PointSize',55,'StdColor','none',...
    'Whiskers','lines','WhiskerLineWidth',2.5,'MarkerEdgeColor',Colors);
box off
hold on
plot(x', y','Color',[0 0 0 0.1],'LineWidth',2)

xlim([0.5 4.5])
ylim([0.02 10])
box off
set(gca,'yscale','log')
ylabel('Firing Rate (Hz)')
set(gca,'YTickLabel',[0.1, 1, 10])
set(gca,'XTick',1:4)
set(gca,'XTickLabel',{'Active','Quiet','NREM','REM'},'XTickLabelRotation',45)
title(['Cell Avg. Fr by State -- All Anims'],'FontSize',16)    


% plot average CV in the 4 states for each cell
figure('Position',[0.2375 0.4500 0.1948 0.4019])

data_toplot = [AW_meanCVbycell,...
    QW_meanCVbycell,...
    NREM_meanCVbycell,...
    REM_meanCVbycell];

Colors = [c_aw;c_qw;c_nrem; c_rem];
[s_obj, x, y] = UnivarScatter_ATP(data_toplot,'BoxType','SEM',...
    'Width',1,'Compression',500,'MarkerFaceColor',Colors,...
    'PointSize',55,'StdColor','none',...
    'Whiskers','lines','WhiskerLineWidth',2.5,'MarkerEdgeColor',Colors);
box off
hold on
plot(x', y','Color',[0 0 0 0.1],'LineWidth',2)

xlim([0.5 4.5])
ylim([0.75 1.75])
box off
ylabel('CV of ISI')
set(gca,'XTick',1:4)
set(gca,'XTickLabel',{'Active','Quiet','NREM','REM'},'XTickLabelRotation',45)
title(['Cell Avg. CV by State -- All Anims - posthunt hrs'],'FontSize',10)    


[p_13,h,stats] = signrank(data_toplot(:,1),data_toplot(:,3));

state_strs = {'Active','Quiet','NREM','REM'};
pairs = nchoosek(1:4,2);
for i = 1:length(pairs)
    pair_i = pairs(i,:);
    [p_i,h,stats] = signrank(data_toplot(:,pair_i(1)),data_toplot(:,pair_i(2)));
    disp([state_strs{pair_i(1)},'-',state_strs{pair_i(2)},' p=',num2str(p_i*length(pairs))])
end


[p,tbl,stats] = anova1(data_toplot)

figure
[c] = multcompare(stats,'display','on')

%% BL v hunt avg
% plot average FR in the 4 states for each cell
figure('Position',[0.0646 0.2611 0.2323 0.5981])

data_toplot = [AW_meanFRbycell_BL,AW_meanFRbycell_hunt,...
    QW_meanFRbycell_BL,QW_meanFRbycell_hunt,...
    NREM_meanFRbycell_BL,NREM_meanFRbycell_hunt,...
    REM_meanFRbycell_BL,REM_meanFRbycell_hunt];

Colors = [c_aw_cont; c_aw; c_qw_cont; c_qw; c_nrem_cont; c_nrem; c_rem_cont; c_rem];
[s_obj, x, y] = UnivarScatter_ATP(data_toplot,'BoxType','SEM',...
    'Width',1,'Compression',500,'MarkerFaceColor',Colors,...
    'PointSize',55,'StdColor','none',...
    'Whiskers','lines','WhiskerLineWidth',2.5,'MarkerEdgeColor',Colors);
box off
hold on
% plot(x', y','Color',[0 0 0 0.1],'LineWidth',2)

for i = 1:2:size(x,2)
    plot(x(:,i:i+1)', y(:,i:i+1)','Color',[0 0 0 0.1],'LineWidth',2)
end


xlim([0.5 8.5])
ylim([0.02 10])
box off
set(gca,'yscale','log')
ylabel('Firing Rate (Hz)')
set(gca,'YTickLabel',[0.1, 1, 10])
set(gca,'XTick',[1.5, 3.5, 5.5, 7.5])
set(gca,'XTickLabel',{'Active','Quiet','NREM','REM'},'XTickLabelRotation',45)
title(['Cell Avg. Fr by State -- All Anims - BLvHunt'],'FontSize',16)    


% plot average CV in the 4 states for each cell
figure('Position',[0.0646 0.2611 0.2323 0.5981])

% data_toplot = [AW_meanCVbycell_BL,...
%     QW_meanCVbycell_BL,...
%     NREM_meanCVbycell_BL,...
%     REM_meanCVbycell_BL];

data_toplot = [AW_meanCVbycell_BL,AW_meanCVbycell_hunt,...
    QW_meanCVbycell_BL,QW_meanCVbycell_hunt,...
    NREM_meanCVbycell_BL,NREM_meanCVbycell_hunt,...
    REM_meanCVbycell_BL,REM_meanCVbycell_hunt];

Colors = [c_aw_cont; c_aw; c_qw_cont; c_qw; c_nrem_cont; c_nrem; c_rem_cont; c_rem];

% Colors = [c_aw;c_qw;c_nrem; c_rem];
[x, y] = UnivarScatter_ATP(data_toplot,'BoxType','SEM',...
    'Width',1,'Compression',500,'MarkerFaceColor',Colors,...
    'PointSize',55,'StdColor','none',...
    'Whiskers','lines','WhiskerLineWidth',2.5,'MarkerEdgeColor',Colors);
box off
hold on
% plot(x', y','Color',[0 0 0 0.1],'LineWidth',2)

for i = 1:2:size(x,2)
    plot(x(:,i:i+1)', y(:,i:i+1)','Color',[0 0 0 0.1],'LineWidth',2)
end

xlim([0.5 8.5])
ylim([0.5 2])
box off
ylabel('CV of ISI')
set(gca,'XTick',[1.5, 3.5, 5.5, 7.5])
% legend({'BL Active','Post hunt Active'})
set(gca,'XTickLabel',{'Active','Quiet','NREM','REM'},'XTickLabelRotation',45)
title(['Cell Avg. CV by State -- All Anims - BLvHunt'],'FontSize',16)    

%% post hunt avg

% plot average FR in the 4 states for each cell
figure('Position',[0.0646 0.2611 0.2323 0.5981])

data_toplot = [AW_meanFRbycell_hunt,...
    QW_meanFRbycell_hunt,...
    NREM_meanFRbycell_hunt,...
    REM_meanFRbycell_hunt];

Colors = [c_aw;c_qw;c_nrem; c_rem];
[s_obj, x, y] = UnivarScatter_ATP(data_toplot,'BoxType','SEM',...
    'Width',1,'Compression',500,'MarkerFaceColor',Colors,...
    'PointSize',55,'StdColor','none',...
    'Whiskers','lines','WhiskerLineWidth',2.5,'MarkerEdgeColor',Colors);
box off
hold on
plot(x', y','Color',[0 0 0 0.1],'LineWidth',2)

xlim([0.5 4.5])
ylim([0.02 10])
box off
set(gca,'yscale','log')
ylabel('Firing Rate (Hz)')
set(gca,'YTickLabel',[0.1, 1, 10])
set(gca,'XTick',1:4)
set(gca,'XTickLabel',{'Active','Quiet','NREM','REM'},'XTickLabelRotation',45)
title(['Cell Avg. Fr by State -- All Anims - hunt'],'FontSize',16)    


% plot average CV in the 4 states for each cell
figure('Position',[0.0646 0.2611 0.2323 0.5981])

data_toplot = [AW_meanCVbycell_hunt,...
    QW_meanCVbycell_hunt,...
    NREM_meanCVbycell_hunt,...
    REM_meanCVbycell_hunt];

Colors = [c_aw;c_qw;c_nrem; c_rem];
[x, y] = UnivarScatter_ATP(data_toplot,'BoxType','SEM',...
    'Width',1,'Compression',500,'MarkerFaceColor',Colors,...
    'PointSize',55,'StdColor','none',...
    'Whiskers','lines','WhiskerLineWidth',2.5,'MarkerEdgeColor',Colors);
box off
hold on
plot(x', y','Color',[0 0 0 0.1],'LineWidth',2)

xlim([0.5 4.5])
ylim([0.5 2])
box off
ylabel('CV of ISI')
set(gca,'XTick',1:4)
set(gca,'XTickLabel',{'Active','Quiet','NREM','REM'},'XTickLabelRotation',45)
title(['Cell Avg. CV by State -- All Anims - hunt'],'FontSize',16)    


%% hunt norm avgs
% plot average FR in the 4 states for each cell
figure('Position',[0.1062 0.5333 0.2031 0.3500])

data_toplot = [AW_meanFRbycell_hunt./AW_meanFRbycell_BL,...
    QW_meanFRbycell_hunt./QW_meanFRbycell_BL,...
    NREM_meanFRbycell_hunt./NREM_meanFRbycell_BL,...
    REM_meanFRbycell_hunt./REM_meanFRbycell_BL];


Colors = [c_aw;c_qw;c_nrem; c_rem];
[s_obj, x, y] = UnivarScatter_ATP(data_toplot,'BoxType','SEM',...
    'Width',1,'Compression',500,'MarkerFaceColor',Colors,...
    'PointSize',55,'StdColor','none',...
    'Whiskers','lines','WhiskerLineWidth',2.5,'MarkerEdgeColor',Colors);
box off
hold on
plot(x', y','Color',[0 0 0 0.1],'LineWidth',2)

xlim([0.5 4.5])
ylim([0.15 35])

r1 = refline(0,1);
r1.LineWidth = 1.5;
r1.LineStyle = '--';
r1.Color = 'k';

box off
set(gca,'yscale','log')
ylabel('Firing Rate (Norm. to Baseline)')
set(gca,'YTickLabel',[1, 10])
set(gca,'XTick',1:4)
set(gca,'XTickLabel',{'Active','Quiet','NREM','REM'},'XTickLabelRotation',45)
title(['Cell Avg. Fr by State -- All Anims - hunt norm - last 12 hrs'],'FontSize',12)    


if strcmp(test_type,'ttest')
    [~,p_aw] = ttest(data_toplot(:,1)-1);
    [~,p_qw] = ttest(data_toplot(:,2)-1);
    [~,p_nrem] = ttest(data_toplot(:,3)-1);
    [~,p_rem] = ttest(data_toplot(:,4)-1);
% 
else
    [p_aw,h,stats] = signrank(data_toplot(:,1)-1);
    [p_qw,h,stats] = signrank(data_toplot(:,2)-1);
    [p_nrem,h,stats] = signrank(data_toplot(:,3)-1);
    [p_rem,h,stats] = signrank(data_toplot(:,4)-1);
end

p_aw = p_aw*4;
p_qw = p_qw*4;
p_nrem = p_nrem*4;
p_rem = p_rem*4;
% [a_aw,fsz_aw] = get_asterisks_from_pval(p_aw);
% [a_qw,fsz_qw] = get_asterisks_from_pval(p_qw);
% [a_nrem,fsz_nrem] = get_asterisks_from_pval(p_nrem);

disp(['p-values  --  AW p=',num2str(p_aw),', QW p=',num2str(p_qw),', NREM p=',num2str(p_nrem),', REM p=',num2str(p_rem)])


% plot average CV in the 4 states for each cell
figure('Position',[0.2375 0.4500 0.1948 0.4019])

data_toplot = [AW_meanCVbycell_hunt./AW_meanCVbycell_BL,...
    QW_meanCVbycell_hunt./QW_meanCVbycell_BL,...
    NREM_meanCVbycell_hunt./NREM_meanCVbycell_BL,...
    REM_meanCVbycell_hunt./REM_meanCVbycell_BL];

Colors = [c_aw;c_qw;c_nrem; c_rem];
[s_obj, x, y] = UnivarScatter_ATP(data_toplot,'BoxType','SEM',...
    'Width',1,'Compression',500,'MarkerFaceColor',Colors,...
    'PointSize',55,'StdColor','none',...
    'Whiskers','lines','WhiskerLineWidth',2.5,'MarkerEdgeColor',Colors);
box off
hold on
plot(x', y','Color',[0 0 0 0.1],'LineWidth',2)

xlim([0.5 4.5])
ylim([0.5 1.8])

r1 = refline(0,1);
r1.LineWidth = 1.5;
r1.LineStyle = '--';
r1.Color = 'k';

box off
ylabel('CV of ISI')
set(gca,'XTick',1:4)
set(gca,'XTickLabel',{'Active','Quiet','NREM','REM'},'XTickLabelRotation',45)
title(['Cell Avg. CV by State -- All Anims - hunt norm - last 12 hrs'],'FontSize',12)    


if strcmp(test_type,'ttest')
    [~,p_aw] = ttest(data_toplot(:,1)-1);
    [~,p_qw] = ttest(data_toplot(:,2)-1);
    [~,p_nrem] = ttest(data_toplot(:,3)-1);
    [~,p_rem] = ttest(data_toplot(:,4)-1);
% 
else
    [p_aw,h,stats] = signrank(data_toplot(:,1)-1);
    [p_qw,h,stats] = signrank(data_toplot(:,2)-1);
    [p_nrem,h,stats] = signrank(data_toplot(:,3)-1);
    [p_rem,h,stats] = signrank(data_toplot(:,4)-1);
end

p_aw = p_aw*4;
p_qw = p_qw*4;
p_nrem = p_nrem*4;
p_rem = p_rem*4;
[a_aw,fsz_aw] = get_asterisks_from_pval(p_aw);
[a_qw,fsz_qw] = get_asterisks_from_pval(p_qw);
[a_nrem,fsz_nrem] = get_asterisks_from_pval(p_nrem);
[a_rem,fsz_nrem] = get_asterisks_from_pval(p_rem);

disp(['p-values  --  AW=',num2str(p_aw),'  QW=',num2str(p_qw),'  NREM=',num2str(p_nrem),'  REM=',num2str(p_rem)])


%%