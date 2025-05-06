
% This script requires: Cell Cluster Structure, behavioral segment
% matrices in the cell structure

% Written by Brian Cary

% cochrans Q test source:
% Jos (10584) (2024). COCHRAN Q TEST (https://www.mathworks.com/matlabcentral/fileexchange/16753-cochran-q-test), MATLAB Central File Exchange. Retrieved December 16, 2024. 


mfilePath = mfilename('fullpath');
if isempty(mfilePath)
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
cd(fileparts(mfilePath))

filepath = "path/to/huntstructure";

load(filepath);

setFigureDefaults
set(0,'defaultFigureUnit','pixels');
set(0,'defaultFigurePosition',[28 530 754 347]);

MASTER = HuntCell.MASTER;

STATES = HuntCell.STATES;
DAYSTART = STATES.DAYSTART;
anim_fields = fieldnames(DAYSTART);
anim_fields(strcmp(anim_fields,'DAYSTART')) = [];
n_anims = numel(anim_fields); % minus daystart field

cell_anims = {MASTER.animal};

already_corrected = 0;
ontime_hr_YES = 1;

% colors for plotting
c_rem   = [25 181 149]./255;
c_nrem  = [131 49 146]./255;
c_aw    = [201 28 101]./255;
c_qw    = [247 148 41]./255;

c_pur = [19 114 186]./255;
c_cons = [255 240 195]./255;
c_both = [106 189 69]./255;

% correct on/off time units to seconds
if ontime_hr_YES && ~already_corrected
    for i = 1:length(MASTER)
        MASTER(i).onTime = MASTER(i).onTime*3600;
        MASTER(i).offTime = MASTER(i).offTime*3600;
    end

    already_corrected = 1;
end


%% Main Loop cell responses to hunt times

Fig_ON_Lines = 0;
Fig_ON_Lines_State = 0;
Fig_ON_Resp = 0;
Run_BOOT_ON = 0;
Fig_ON_Boot = 0;
Fig_ON_TC_cell = 0;
Fig_ON_TC_anim = 0;

Fig_Boot_ManHandles = 0;
man_boot_y_max = 2.5;

celltype = 'RSU'; %RSU or pFS

save_fig_ON = 0;
save_path = 'C:\Users\bcary\Documents\Work\LabCode\DanCode\figs\';
save_types = {'.fig','.pdf','.tiff'};

smooth_ON = 1; % on for timecourses
fr_bin_sec = 0.1; % default=0.1
gauss_wind_sec = 3; % default=3 for RSU, 5 for FS
win_size = round(gauss_wind_sec/fr_bin_sec); % make x sec gauss wind
gauss_wind = gausswin(win_size);
gauss_wind = gauss_wind/sum(gauss_wind);

fr_change_calc = 2; %1 for fractional (purFR/iciFR), 2 for % change (pur-ici)/ici;

fr_bin_sec_tc = 0.2; % default=0.2
gauss_wind_sec = 3;
win_size = round(gauss_wind_sec/fr_bin_sec_tc); % make x sec gauss wind
gauss_wind_tc = gausswin(win_size);
gauss_wind_tc = gauss_wind_tc/sum(gauss_wind_tc);

min_behFr_ici = 0.03; %0.05 for FS 0.03 for RSU

boot_ext_hrs = 1;
% boot_hr_range = 4;
num_boots = 500;
boot_low_prc = 5;
boot_upp_prc = 95;
kick_out_ON = 1; % removes ICI FRs if they are below a certain low bootstrap percentile or min thresh
min_boot_PurFR = 0.05;

center_beh = 2;
tc_pre_sec = 15;
tc_post_sec = 15;
tc_baseline_ICI = 0;
tc_baseline_sec = 8;
tc_min_base_prop = 0.5;
tc_min_base_fr = 0.1;
tc_time = -tc_pre_sec:fr_bin_sec_tc:(tc_post_sec - fr_bin_sec_tc);

if Run_BOOT_ON == 1
    full_boot_sigs = {};
    full_boot_percs = {};
end
full_cap_timecourses = {};
full_behFRs = {};
for anim_i = 1:n_anims
    
    %% init
%     anim_states = STATES.(anim_fields{anim_i});
    anim_daystart = DAYSTART.(anim_fields{anim_i});
%     state_start = anim_states(1,2);

    anim_inds = strcmp(cell_anims,anim_fields{anim_i});
    anim_cells = MASTER(anim_inds);

%     if isempty(anim_cells)
%         continue
%     end

    expt_start = anim_cells(1).EXPTSTART;
    trem = anim_cells(1).trem;

    
    % format beh segs
    Beh_Seg_str = {};
    for seg = 1:5
        try
            Beh_Seg_str{seg} = anim_cells(1).(['behavioralsegments',num2str(seg)]);
        end
    end

    %% Load in behavioral start stop times into structure
    ICI_OnOff = {};
    Pur_OnOff = {};
    Cons_OnOff = {};

    for beh_num = 1:length(Beh_Seg_str)
    
        beh_seg = Beh_Seg_str{beh_num};
            
        if ~isempty(beh_seg)
    
            ICI_idx = [];
            Pur_idx = [];
            Cons_idx = [];
            
            ICI_idx = find(beh_seg(:,3) == 1);
            ICI_idx = [find(beh_seg(:,3) == 5); ICI_idx];
            ICI_idx = sort(ICI_idx);
            ICI_start = beh_seg(ICI_idx,1);
            ICI_stop = beh_seg(ICI_idx,2);
        
            Pur_idx = find(beh_seg(:,3) == 3);
            Pur_start = nan([1,length(ICI_start)])';
            Pur_start(1:length(Pur_idx)) = beh_seg(Pur_idx,1);
            Pur_stop = nan([1,length(ICI_start)])';
            Pur_stop(1:length(Pur_idx)) = beh_seg(Pur_idx,2);
        
            Cons_idx = find(beh_seg(:,3) == 4);
            Cons_start = nan([1,length(ICI_start)])';
            Cons_start(1:length(Cons_idx)) = beh_seg(Cons_idx,1);
            Cons_stop = nan([1,length(ICI_start)])';
            Cons_stop(1:length(Cons_idx)) = beh_seg(Cons_idx,2);
    
            ICI_durs = ICI_stop - ICI_start;
            Pur_durs = Pur_stop - Pur_start;
            Cons_durs = Cons_stop - Cons_start;
            
            disp(['Behavioral Session ',num2str(beh_num)])
            disp(['Avg ICI dur=',num2str(nanmean(ICI_durs)),'   std=',num2str(nanstd(ICI_durs))])
            disp(['Avg Pur dur=',num2str(nanmean(Pur_durs)),'   std=',num2str(nanstd(Pur_durs))])
            disp(['Avg Cons. dur=',num2str(nanmean(Cons_durs)),'   std=',num2str(nanstd(Cons_durs))])
            
            ICI_OnOff{beh_num,1} = ICI_start;
            ICI_OnOff{beh_num,2} = ICI_stop;
            Pur_OnOff{beh_num,1} = Pur_start;
            Pur_OnOff{beh_num,2} = Pur_stop;    
            Cons_OnOff{beh_num,1} = Cons_start;
            Cons_OnOff{beh_num,2} = Cons_stop;   
            
            if length(ICI_start) ~= length(Pur_start)
                disp('different number of ICI and Pur starts')
            end
            if length(ICI_stop) ~= length(Pur_stop)
                disp('different number of ICI and Pur stops')
            end  
        end
    end
        
    %% cell loop
    anim_boot_sigs = {};
    anim_boot_percs = {};
    anim_cap_timecourses = {};
    anim_behFrs = {};
    anim_boot_sigs(1:size(ICI_OnOff,1)) = {nan([length(anim_cells),2])};
    anim_boot_percs(1:size(ICI_OnOff,1)) = {nan([length(anim_cells),2])};
    anim_behFrs(1:size(ICI_OnOff,1)) = {nan([length(anim_cells),2])};
    anim_cap_timecourses(1:size(ICI_OnOff,1)) = {nan([length(anim_cells),length(tc_time)])};
    for cell_i = 1:length(anim_cells)
        
        if strcmp(anim_cells(cell_i).CellType,celltype) ~= 1
            continue
        end

        disp(['Processing cell: ',num2str(cell_i), '  anim num=',num2str(anim_i)])
        tic
        
        %% initialize 
        all_sp_times = anim_cells(cell_i).time;
        onTimes = anim_cells(cell_i).onTime;
        offTimes = anim_cells(cell_i).offTime;
            
        cell_sp_times = [];
        num_ontimes = length(onTimes);
        for ii = 1:num_ontimes
            add_sps = all_sp_times(all_sp_times > onTimes(ii)...
                & all_sp_times < offTimes(ii));
            cell_sp_times = [cell_sp_times; add_sps];
        end
    
        
        if Fig_ON_Lines == 1
            fig_lines1 = figure('Position',[50 115 1142 809]); hold on;
        end

        if Fig_ON_Lines_State == 1
            fig_lines2 = figure('Position',[75 115 1142 809]); hold on;
        end
        
        if Fig_ON_Resp == 1
            fig_resp1 = figure('Position',[44 610 1113 370]);
        end
        
        if Fig_ON_Boot == 1
            if Fig_Boot_ManHandles ~= 1
                fig_resp2 = figure('Position',[143 655 908 315]);
            end
        end

        if Fig_ON_TC_cell == 1
            fig_tc1 = figure('Position',[65 665 809 315]);
        end

        %% behavior loop
        for beh_seg = 1:size(ICI_OnOff,1)
            
            %% init
            % get start of diff beh's in seconds from structures
            ICI_starts = ICI_OnOff{beh_seg,1};
            ICI_ends = ICI_OnOff{beh_seg,2};
            Pur_starts = Pur_OnOff{beh_seg,1};
            Pur_ends = Pur_OnOff{beh_seg,2};
            Cons_starts = Cons_OnOff{beh_seg,1};
            Cons_ends = Cons_OnOff{beh_seg,2};        
            beh_end = max([ICI_ends;Pur_ends;Cons_ends]);

            % extra check for beh inds
            for crick_i = 1:length(ICI_starts)
                if ICI_ends(crick_i) < ICI_starts(crick_i)
                    disp(['Bad ICI_end label: ','crick=',num2str(crick_i)])
                    ICI_ends(crick_i) = Pur_starts(crick_i) - 1;
                end

                if Pur_ends(crick_i) < Pur_starts(crick_i)
                    disp(['Bad Pur_end label: ','crick=',num2str(crick_i)])
                    Pur_ends(crick_i) = Cons_starts(crick_i);
                end

                if Cons_ends(crick_i) < Cons_starts(crick_i)
                    disp(['Bad Cons_ends label: ','crick=',num2str(crick_i)])
                    if crick_i ~= length(ICI_starts)
                        Cons_ends(crick_i) = ICI_starts(crick_i+1);
                    end
                end                
            end

            % use histcounts to create FR histogram over time
            edges = ICI_starts(1):fr_bin_sec:beh_end;
            [N,~] = histcounts(cell_sp_times,edges);
            cell_fr = N/fr_bin_sec;
            
            % apply gauss filt to cell FR
            if smooth_ON == 1
                cell_fr = filter(gauss_wind,1,cell_fr);
            end
                
            % re-center time bins
            time_bin = edges(2:end) - fr_bin_sec/2;

%             % create state vec
%             %1 = REM, 2 = NREM, 4 = AW, 5 = QW
%             state_times = anim_states(:,2);
%             state_times = state_times - expt_start + trem;
%             state_vec_resp = nan(1,length(time_bin));
%             for s = 1:length(state_times)-1
%                 s_inds = time_bin > state_times(s) & time_bin <= state_times(s+1);
%                 state_vec_resp(s_inds) = anim_states(s,1);
%             end
%             state_vec_resp = state_vec_resp(2:end);    

            %% Calculate avg FR responses to diff behavioral segments
            
            ICI_FRs = [];
            Pur_FRs = [];
            Cons_FRs = [];
            for crick = 1:length(ICI_starts)
                ICI_FRs(crick) = nanmean(cell_fr(time_bin > ICI_starts(crick) & time_bin < ICI_ends(crick)));
                Pur_FRs(crick) = nanmean(cell_fr(time_bin >= Pur_starts(crick) & time_bin < Pur_ends(crick)));
                Cons_FRs(crick) = nanmean(cell_fr(time_bin >= Cons_starts(crick) & time_bin < Cons_ends(crick)));
            end
            
            ICI_FRs(ICI_FRs < min_behFr_ici) = NaN;
            
            if fr_change_calc == 1
                % normalized to each prior ICI
                anim_behFrs{beh_seg}(cell_i,:) = [nanmean(Pur_FRs./ICI_FRs), nanmean(Cons_FRs./ICI_FRs)];
    
                % normalizing the averages to one another
    %             anim_behFrs{beh_seg}(cell_i,:) = [nanmean(Pur_FRs)./nanmean(ICI_FRs), nanmean(Cons_FRs)./nanmean(ICI_FRs)];
            elseif fr_change_calc ==2
                anim_behFrs{beh_seg}(cell_i,:) = [100*nanmean((Pur_FRs - ICI_FRs)./ICI_FRs),...
                                                  100*nanmean((Cons_FRs - ICI_FRs)./ICI_FRs)];

            end

            %% Plotting lines
            if Fig_ON_Lines == 1
                
                % annoying way I have to find indices for differnt behaviors to
                % plot them as different colors in the same line plot
                ICI_idx = zeros(1,length(time_bin));
                Pur_idx = zeros(1,length(time_bin));
                Cons_idx = zeros(1,length(time_bin));
                for crick = 1:length(ICI_starts)
                    new_idx = time_bin > ICI_starts(crick) & time_bin < ICI_ends(crick);
                    ICI_idx = ICI_idx + new_idx;
    
                    new_idx = time_bin >= Pur_starts(crick) & time_bin < Pur_ends(crick);
                    Pur_idx = Pur_idx + new_idx;      
    
                    new_idx = time_bin >= Cons_starts(crick) & time_bin < Cons_ends(crick);
                    Cons_idx = Cons_idx + new_idx;   
                end
    
                ax = subplot(size(ICI_OnOff,1),1,beh_seg,'Parent',fig_lines1);
%                 
%                 figure('Position',[9 686 565 291])
%                 ax = gca;
                plot(ax,time_bin,cell_fr,'k');
                hold(ax,'on')       
                ylabel(ax,'Firing Rate (Hz)')
                xlabel(ax,'Time (s)')
                title_str = ['Cell: ',num2str(cell_i),'  Hunt Sesh: ',num2str(beh_seg),...
                    '  ICI FR=',num2str(nanmean(ICI_FRs)),'  PurNormFR=',num2str(nanmean(Pur_FRs./ICI_FRs))];
                title(ax,title_str,'FontSize',10)        
    
                plot(ax,Pur_starts, max(cell_fr)*ones(1,length(Pur_starts)),'g*')    
                plot(ax,Pur_ends, max(cell_fr)*ones(1,length(Pur_ends)),'r*')  
    
                ICI_fr = cell_fr;
                ICI_fr(~ICI_idx) = NaN;
                Pur_fr = cell_fr;
                Pur_fr(~Pur_idx) = NaN;
                Cons_fr = cell_fr;
                Cons_fr(~Cons_idx) = NaN;
                
                plot(ax,time_bin,ICI_fr,'b');
                plot(ax,time_bin,Pur_fr,'g');
                plot(ax,time_bin,Cons_fr,'m');
                box(ax,'off')
    
    %             xlim([277970, 278420])
    %             xlim([274530, 274530+450])
    %             ylim([0, 20])
            end
            
            
            if Fig_ON_Lines_State == 1
                
                ax = subplot(size(ICI_OnOff,1),1,beh_seg,'Parent',fig_lines2);
                
    %             figure('Position',[9 686 565 291])
    %             ax = gca;
                plot(ax,time_bin,cell_fr,'k');
                hold(ax,'on')       
                ylabel(ax,'Firing Rate (Hz)')
                xlabel(ax,'Time (s)')
                title_str = ['Cell: ',num2str(cell_i),'  Hunt Sesh: ',num2str(beh_seg)];
                title(ax,title_str)        
    
                plot(ax,Pur_starts, max(cell_fr)*ones(1,length(Pur_starts)),'g*')    
                plot(ax,Pur_ends, max(cell_fr)*ones(1,length(Pur_ends)),'r*')  
    
                nrem_fr = cell_fr;
                rem_fr = cell_fr;
                qw_fr = cell_fr;
                aw_fr = cell_fr;
                nrem_fr(state_vec_resp ~= 2) = NaN;
                rem_fr(state_vec_resp ~= 1) = NaN;
                qw_fr(state_vec_resp ~= 5) = NaN;
                aw_fr(state_vec_resp ~= 4) = NaN;
                
                line_objs1= plot(ax,time_bin,nrem_fr,'color',c_nrem);
                line_objs2 = plot(ax,time_bin,rem_fr,'color',c_rem);  
                line_objs3 = plot(ax,time_bin,qw_fr,'color',c_qw);
                line_objs4 = plot(ax,time_bin,aw_fr,'color',c_aw);
    
                title_str = ['Cell: ',num2str(cell_i)];
                title(ax,title_str)        
                    
                box(ax,'off')
                leg_names = {'nrem','rem','qw','aw'};
                legend(ax,[line_objs1,line_objs2,line_objs3,line_objs4],leg_names);
                box(ax,'off')
    
    %             xlim([277970, 278420])
    %             xlim([274530, 274530+450])
    %             ylim([0, 20])
            end
            

            %% Bootstrap 
            if Run_BOOT_ON == 1 && nanmean(Pur_FRs) > min_boot_PurFR
                tic
                %% calc bootstrap
                % use histcounts to create FR histogram over time
                boot_ext_inds = ceil(boot_ext_hrs*3600/fr_bin_sec);
                beh_st = ICI_starts(1);
                beh_end = max([ICI_ends;Pur_ends;Cons_ends]);
                beh_dur = (beh_end - beh_st);
%                 beh_mid = ICI_starts(1) + beh_dur/2;
%                 boot_st = beh_mid-boot_hr_range*3600/2;
%                 edges = boot_st:fr_bin_sec:(beh_mid+boot_hr_range*3600/2);

                boot_st = beh_st - boot_ext_hrs*3600;
                boot_end = beh_end + boot_ext_hrs*3600;
                edges = boot_st:fr_bin_sec:boot_end;
                [N,~] = histcounts(cell_sp_times,edges);
                cell_fr_boot = N/fr_bin_sec;        
    
                % apply gauss filt to cell FR
                if smooth_ON == 1
                    cell_fr_boot = filter(gauss_wind,1,cell_fr_boot);
                end        
                time_bin = round(edges(2:end) - fr_bin_sec/2, 2);
    
%                 % create state vec
%                 %1 = REM, 2 = NREM, 4 = AW, 5 = QW
%                 state_times = anim_states(:,2);
%                 state_times = state_times - expt_start + trem;
%                 state_vec_boot = nan(1,length(time_bin));
%                 for s = 1:length(state_times)-1
%                     s_inds = time_bin > state_times(s) & time_bin <= state_times(s+1);
%                     state_vec_boot(s_inds) = anim_states(s,1);
%                 end
%                 state_vec_boot = state_vec_boot(2:end);    

                cell_fr_mat = repmat(cell_fr_boot',[1 num_boots]);
    
                ICI_st_idx = round((ICI_starts-beh_st)*(1/fr_bin_sec))+1;
                ICI_end_idx = round((ICI_ends-beh_st)*(1/fr_bin_sec))+1;
                Pur_st_idx = round((Pur_starts-beh_st)*(1/fr_bin_sec))+1;
                Pur_end_idx = round((Pur_ends-beh_st)*(1/fr_bin_sec))+1;
                Cons_st_idx = round((Cons_starts-beh_st)*(1/fr_bin_sec))+1;
                Cons_end_idx = round((Cons_ends-beh_st)*(1/fr_bin_sec))+1;

%                 ICI_st_idx = round((ICI_starts-boot_st)*(1/fr_bin_sec));
%                 ICI_end_idx = round((ICI_ends-boot_st)*(1/fr_bin_sec));
%                 Pur_st_idx = round((Pur_starts-boot_st)*(1/fr_bin_sec));
%                 Pur_end_idx = round((Pur_ends-boot_st)*(1/fr_bin_sec));
%                 Cons_st_idx = round((Cons_starts-boot_st)*(1/fr_bin_sec));
%                 Cons_end_idx = round((Cons_ends-boot_st)*(1/fr_bin_sec));        
    
%                 sample_time = time_bin(1:end-round(beh_dur*(1/fr_bin_sec)));
%                 shifts = randsample(1:length(sample_time),num_boots) - ICI_st_idx(1);

                shifts = randsample(1:(boot_ext_inds*2-2),num_boots);

                ICI_FRs_boot = zeros(length(ICI_starts),num_boots);
                Pur_FRs_boot = zeros(length(ICI_starts),num_boots);
                Cons_FRs_boot = zeros(length(ICI_starts),num_boots);
                for crick = 1:length(ICI_starts)

                    if ~isnan(ICI_st_idx(crick))
                        ICI_boot_idx = zeros(length(cell_fr_boot),num_boots);
                        Pur_boot_idx = zeros(length(cell_fr_boot),num_boots);
                        Cons_boot_idx = zeros(length(cell_fr_boot),num_boots);
        
                        for ii = 1:num_boots
                            ici_inds = (ICI_st_idx(crick) + shifts(ii)):(ICI_end_idx(crick) + shifts(ii));
                            ICI_boot_idx(ici_inds,ii) = 1;
                            
                            if ~isnan(Pur_st_idx(crick))
                                pur_inds = (Pur_st_idx(crick) + shifts(ii)):(Pur_end_idx(crick) + shifts(ii));
                                Pur_boot_idx(pur_inds,ii) = 1;   
                                
                                cons_inds = (Cons_st_idx(crick) + shifts(ii)):(Cons_end_idx(crick) + shifts(ii));
                                Cons_boot_idx(cons_inds,ii) = 1;   
                            end
                        end
        
                        ICI_boot_frs = reshape(cell_fr_mat(logical(ICI_boot_idx)),[length(ici_inds) num_boots]);
                        ICI_FRs_boot(crick,:) = nanmean(ICI_boot_frs,1);
                        if ~isnan(Pur_st_idx(crick))
                            Pur_boot_frs = reshape(cell_fr_mat(logical(Pur_boot_idx)),[length(pur_inds) num_boots]);
                            Pur_FRs_boot(crick,:) = nanmean(Pur_boot_frs,1);            
                            Cons_boot_frs = reshape(cell_fr_mat(logical(Cons_boot_idx)),[length(cons_inds) num_boots]);
                            Cons_FRs_boot(crick,:) = nanmean(Cons_boot_frs,1); 
                        end
                    else
                        ICI_FRs_boot(crick,:) = nan([1,num_boots]);
                        Pur_FRs_boot(crick,:) = nan([1,num_boots]);
                        Cons_FRs_boot(crick,:) = nan([1,num_boots]);
                    end
                end 
    
                ici_boot_mean = nanmean(ICI_FRs_boot,1);
                pur_boot_mean = nanmean(Pur_FRs_boot,1);
                cons_boot_mean = nanmean(Cons_FRs_boot,1);
                
                if Fig_ON_Boot == 1
                    if Fig_Boot_ManHandles~=1
                        % plot results
                        ax_boot = subplot(1,size(ICI_OnOff,1),beh_seg,'Parent',fig_resp2);
                        hold(ax_boot,'on')

                    else
    
                        %manual boot plotting
                        fig_resp2 = figure('Position',[252 267 320 385]);
                        ax_boot = axes('parent',fig_resp2);
                        hold(ax_boot,'on')
                    end
                end

                if kick_out_ON == 1
%                     out_inds = ICI_FRs < prctile(ici_boot_mean,0.5);
                    out_inds = ICI_FRs < min_boot_PurFR;

                    ICI_FRs(out_inds) = NaN;
                    Pur_FRs(out_inds) = NaN;
                    Cons_FRs(out_inds) = NaN;

                    if sum(out_inds) > 0
                        disp('Kicking out cricks from boot analysis: ')
                        disp(find(out_inds))
                    end
                end
                
%                 ici_boot_n = ici_boot_mean./nanmean(ici_boot_mean);
%                 pur_boot_n = pur_boot_mean./nanmean(ici_boot_mean);
%                 cons_boot_n = cons_boot_mean./nanmean(ici_boot_mean);

                ici_boot_n = ici_boot_mean./(ici_boot_mean);
                pur_boot_n = pur_boot_mean./(ici_boot_mean);
                cons_boot_n = cons_boot_mean./(ici_boot_mean);
                
                p_xs(1,:) = [0.5 1.5 2.5]; 
                p_xs(2,:) = p_xs(1,:) + 0.5;
                p_xs(3,:) = p_xs(2,:);
                p_xs(4,:) = p_xs(1,:);
                p_ys(1:2,:) = [prctile(ici_boot_n,boot_low_prc),...
                    prctile(pur_boot_n,boot_low_prc),...
                    prctile(cons_boot_n,boot_low_prc);...
                    prctile(ici_boot_n,boot_low_prc),...
                    prctile(pur_boot_n,boot_low_prc),...
                    prctile(cons_boot_n,boot_low_prc)];
                p_ys(3:4,:) = [prctile(ici_boot_n,boot_upp_prc),...
                    prctile(pur_boot_n,boot_upp_prc),...
                    prctile(cons_boot_n,boot_upp_prc);...
                    prctile(ici_boot_n,boot_upp_prc),...
                    prctile(pur_boot_n,boot_upp_prc),...
                    prctile(cons_boot_n,boot_upp_prc)];  
                
                x_pts = [ones(1,length(ICI_FRs)); 2*ones(1,length(Pur_FRs)); 3*ones(1,length(Cons_FRs))];
                y_pts = [ICI_FRs./nanmean(ICI_FRs);...
                    Pur_FRs./ICI_FRs;...
                    Cons_FRs./ICI_FRs];   
                
                %% fig boot on
                if Fig_ON_Boot == 1
                    ylabel(ax_boot,'Firing Rate (Norm.)')
                    title_str = {['Anim: ',anim_fields{anim_i}],...
                        ['Cell: ',num2str(cell_i),'  Hunt: ',num2str(beh_seg)]};
                    title(ax_boot,title_str,'FontSize',10);   

                    if Fig_Boot_ManHandles~=1
                        p=patch(ax_boot,p_xs, p_ys,'k'); hold on;
                        set(p,'FaceAlpha',0.15);
                        set(p,'LineStyle','none');   
                    else
                        %manual ylim
                        ylim(ax_boot,[0, man_boot_y_max])

    %                     create way of having embedded sideways histograms of
    %                     boostrap for visualization...
                        xstart=0.45;
                        xend=0.52;
                        ystart=0.1115;
                        yend=0.9;
                        pur_hist_ax = axes('parent',fig_resp2,'position',[xstart ystart xend-xstart yend-ystart]);                 
                        bins = 0:.025:man_boot_y_max;
                        h_pur = histogram(pur_boot_n,bins,'EdgeColor','None','FaceColor',[0 0 0],'FaceAlpha',0.3,...
                            'orientation','horizontal');                    
                        box(pur_hist_ax,'off')
                        set(pur_hist_ax, 'XTick', [], 'XTickLabel', []);
                        set(pur_hist_ax, 'YTick', [], 'YTickLabel', []);
                        set(get(pur_hist_ax, 'XAxis'), 'Visible', 'off');
                        set(get(pur_hist_ax, 'YAxis'), 'Visible', 'off');
                        pur_hist_ax.Color = 'None';
                        ylim(pur_hist_ax,[0, man_boot_y_max])
    
                        xstart=0.68;
                        xend=0.75;
                        ystart=0.1115;
                        yend=0.9;
                        cons_hist_ax = axes('parent',fig_resp2,'position',[xstart ystart xend-xstart yend-ystart]);                 
                        h_cons = histogram(cons_boot_n,bins,'EdgeColor','None','FaceColor',[0 0 0],'FaceAlpha',0.3,...
                            'orientation','horizontal');                    
                        box(cons_hist_ax,'off')
                        set(cons_hist_ax, 'XTick', [], 'XTickLabel', []);
                        set(cons_hist_ax, 'YTick', [], 'YTickLabel', []);
                        set(get(cons_hist_ax, 'XAxis'), 'Visible', 'off');
                        set(get(cons_hist_ax, 'YAxis'), 'Visible', 'off');
                        cons_hist_ax.Color = 'None';
                        ylim(cons_hist_ax,[0, man_boot_y_max])

                        xmax = max([h_pur.Values h_cons.Values]);
                        xlim(cons_hist_ax,[0, xmax])
                        xlim(pur_hist_ax,[0, xmax])
                        
                    end
%                     fig_resp2.Children=flip(fig_resp2.Children);

                    plot(ax_boot,x_pts,y_pts,'Color',[0 0 0 0.1]); hold on;
                                            
                    plot(ax_boot,ones(1,length(ICI_FRs)),y_pts(1,:),'.','Markersize',20,...
                        'Color',[0.5 0.5 0.5]);
                    plot(ax_boot,2*ones(1,length(Pur_FRs)),y_pts(2,:),'.','Markersize',20,...
                        'Color',c_pur);
                    plot(ax_boot,3*ones(1,length(Cons_FRs)),y_pts(3,:),'.','Markersize',20,...
                        'Color',c_cons.*0.88);           
                    
                    plot(ax_boot,[1 2 3],nanmean(y_pts,2),'k','LineWidth',2.5)
                    
                    xticks(ax_boot,[1 2 3])
                    set(ax_boot,'XTickLabel',{'ICI','Pur','Cons'},'XTickLabelRotation',45,'FontSize',10);
                    box(ax_boot, 'off')
                    xlim(ax_boot,[0.5 3.5])
                    
                    maxval = max(y_pts(:));
                    if maxval ~= 0 && ~isnan(maxval)
                        ylim(ax_boot,[0 1.2*maxval])
                    end

                    if Fig_Boot_ManHandles==1
                        ylim(ax_boot,[0, man_boot_y_max])
                    end

                end
                
                %% rest of fig analysis run

                disp(['Pursuit mean: ',num2str(nanmean(y_pts(2,:))),' (boot percentiles: ',...
                    num2str(prctile(pur_boot_n,boot_low_prc)),'-',num2str(prctile(pur_boot_n,boot_upp_prc)),')'])
                disp(['Cons mean: ',num2str(nanmean(y_pts(3,:))),' (boot percentiles: ',...
                    num2str(prctile(cons_boot_n,boot_low_prc)),'-',num2str(prctile(cons_boot_n,boot_upp_prc)),')'])  

                anim_boot_sigs{beh_seg}(cell_i,1) = 0;
                anim_boot_sigs{beh_seg}(cell_i,2) = 0;
                if (nanmean(y_pts(2,:)) > prctile(pur_boot_n,boot_upp_prc))
                    if Fig_ON_Boot == 1
                        text(ax_boot,2,maxval*1.1,'*','Fontsize',20);
                    end
                    anim_boot_sigs{beh_seg}(cell_i,1) = 1;

                elseif (nanmean(y_pts(2,:)) < prctile(pur_boot_n,boot_low_prc))

                    if Fig_ON_Boot == 1
                        text(ax_boot,2,maxval*1.1,'*','Fontsize',20);
                    end

                    anim_boot_sigs{beh_seg}(cell_i,1) = -1;
                end
                if (nanmean(y_pts(3,:)) > prctile(cons_boot_n,boot_upp_prc))
                    if Fig_ON_Boot == 1
                        text(ax_boot,3,maxval*1.1,'*','Fontsize',20);
                    end

                    anim_boot_sigs{beh_seg}(cell_i,2) = 1;

                elseif (nanmean(y_pts(3,:)) < prctile(cons_boot_n,boot_low_prc))

                    if Fig_ON_Boot == 1
                        text(ax_boot,3,maxval*1.1,'*','Fontsize',20);
                    end

                    anim_boot_sigs{beh_seg}(cell_i,2) = -1;
                end            
                    
                if Fig_ON_Boot == 1
                    ylabel(ax_boot,'Firing Rate (Norm.)')
                    title_str = {['Anim: ',anim_fields{anim_i}],...
                        ['Cell: ',num2str(cell_i),'  Hunt: ',num2str(beh_seg)]};
                    title(ax_boot,title_str,'FontSize',10);   
                end
                    
                purs_perc_ind = find(sort(pur_boot_n) > nanmean(Pur_FRs./ICI_FRs), 1, 'first');
                cons_perc_ind = find(sort(cons_boot_n) > nanmean(Cons_FRs./ICI_FRs), 1, 'first');
                
                if isempty(purs_perc_ind)
                    purs_perc_ind = length(pur_boot_n);
                end
                if isempty(cons_perc_ind)
                    cons_perc_ind = length(cons_boot_n);
                end

                purs_perc = 100*purs_perc_ind/length(pur_boot_n);
                cons_perc = 100*cons_perc_ind/length(cons_boot_n);

                anim_boot_percs{beh_seg}(cell_i,:) = [purs_perc, cons_perc];
                
                %%
            end
            drawnow
            toc
            
            %% response plots
            if Fig_ON_Resp == 1
                
                ax_resp = subplot(1,size(ICI_OnOff,1),beh_seg,'Parent',fig_resp1);
                
                x_pts = [ones(1,length(ICI_FRs)); 2*ones(1,length(Pur_FRs)); 3*ones(1,length(Cons_FRs))];
                y_pts = [ICI_FRs; Pur_FRs; Cons_FRs];
                plot(ax_resp,x_pts,y_pts,'Color',[0 0 0 0.15]);
                hold(ax_resp,'on')
                
                plot(ax_resp,[1 2 3],nanmean(y_pts,2),'k','LineWidth',2.5)
                
                plot(ax_resp,ones(1,length(ICI_FRs)),ICI_FRs,'bo');
                plot(ax_resp,2*ones(1,length(Pur_FRs)),Pur_FRs,'go');
                plot(ax_resp,3*ones(1,length(Cons_FRs)),Cons_FRs,'mo');            
                
                xticks(ax_resp,1:3)
                set(ax_resp,'XTickLabel',{'ICI','Pur','Cons'},'XTickLabelRotation',45,'FontSize',10);
                box(ax_resp, 'off')
                xlim(ax_resp,[0.5 3.5])
                
                maxval = max([ICI_FRs Pur_FRs Cons_FRs]);
                if maxval ~= 0
                    ylim(ax_resp,[0 1.2*maxval])
                end
                
                ylabel(ax_resp,'Firing Rate (Hz)')
                title_str = ['Cell: ',num2str(cell_i),'  Hunt Sesh: ',num2str(beh_seg)];
                title(ax_resp,title_str);                  
                
            end
    
            drawnow
        
            %% PSTH timecourses
            
            % use histcounts to create FR histogram over time
            edges_tc = (ICI_starts(1)-tc_pre_sec):fr_bin_sec_tc:(beh_end+tc_post_sec);
            [N,~] = histcounts(cell_sp_times,edges_tc);
            cell_fr_tc = N/fr_bin_sec_tc;
            
            % apply gauss filt to cell FR
            if smooth_ON == 1
                cell_fr_tc = filter(gauss_wind_tc,1,cell_fr_tc);
            end
                
            % re-center time bins
            time_bin_tc = edges_tc(2:end) - fr_bin_sec_tc/2;
            
            if center_beh == 1
                center_times = Pur_starts;
            elseif center_beh == 2
                center_times = Pur_ends;
            end

            cap_FR_timecourse = [];
            cap_FR_timecourse_norm = [];
            for crick_i = 1:length(center_times)
                if ~isnan(center_times(crick_i))
                    start_sec = center_times(crick_i) - tc_pre_sec;
                    end_sec = center_times(crick_i) + tc_post_sec;
    
                    cap_inds = time_bin_tc >= start_sec & time_bin_tc < end_sec;
                    cap_FR_timecourse(crick_i,:) = cell_fr_tc(cap_inds);
    
                    baseline_inds = time_bin_tc >= start_sec & time_bin_tc < (start_sec + tc_baseline_sec);

                    if tc_baseline_ICI == 1
                        base_fr = ICI_FRs(crick_i);
                    else
                        base_fr = nanmean(cell_fr_tc(baseline_inds));
                    end

                    base_prop = sum(baseline_inds)/(tc_baseline_sec/fr_bin_sec_tc);
                    if base_fr > tc_min_base_fr && base_prop > tc_min_base_prop
                        cap_FR_timecourse_norm(crick_i,:) = cell_fr_tc(cap_inds) ./ base_fr;
                    else
                        cap_FR_timecourse_norm(crick_i,:) = nan([1,length(tc_time)]);
                    end
    
                    if ~isempty(find(cap_FR_timecourse_norm==Inf))
                        disp('inf in capture tc')
                        keyboard
                    end
                end

            end
            
            anim_cap_timecourses{beh_seg}(cell_i,:) = nanmean(cap_FR_timecourse_norm,1);

            if Fig_ON_TC_cell == 1
                
                ax_tc = subplot(1,size(ICI_OnOff,1),beh_seg,'Parent',fig_tc1);
                
                tc_data_toplot = cap_FR_timecourse_norm;
                hold(ax_tc,'on')
                alpha = 0.3;
                sh_color = 'g';
                s1 = stdshade(tc_data_toplot,alpha,sh_color,tc_time,1,'sem');   

                maxval = max(nanmean(tc_data_toplot,1));
                if maxval ~= 0
                    ylim(ax_tc,[0 1.2*maxval])
                end

%                 plot([0,0],[0.5, 2],'k')
%                 xlim(ax_tc,[0.5 3.5])
%                 ylim([0.75, 1.5])

                plot([0,0],[0, maxval],'k')
                xlim(ax_tc,[-tc_pre_sec, tc_post_sec])
                
%                 ylabel('Norm FR')
                xlabel('Seconds since capture')                
                ylabel(ax_tc,'Firing Rate (Hz)')
                title_str = ['Cell: ',num2str(cell_i),'  Hunt Sesh: ',num2str(beh_seg)];
                title(ax_tc,title_str);    
                
            end
        
            drawnow

        end
        
        toc
        
        %% save figure
        if save_fig_ON == 1
            
            figtitle = ['Cell',num2str(cell_i),' ',anim_cells(cell_i).CellType];
            try
                for yy = 1:length(save_types)
                    file_savepath = [save_path,filesep,figtitle,save_types{yy}];
                    saveas(fig_resp2,file_savepath)
                end
            catch
                disp('couldnt save -- probably something wrong with title')
                keyboard
            end
    
        end
    
    end

    
    %% anim fig and data store
    if Fig_ON_TC_anim == 1
        figure;
        for beh_i = 1:length(anim_cap_timecourses)
            ax_tc = subplot(1,length(anim_cap_timecourses),beh_i);
            
            hold(ax_tc,'on')
            alpha = 0.3;
            sh_color = 'g';
            s1 = stdshade(anim_cap_timecourses{beh_i},alpha,sh_color,tc_time,1,'sem');    
    
            maxval = max(nanmean(anim_cap_timecourses{beh_i}));
            if maxval ~= 0 && ~isnan(maxval)
                ylim(ax_tc,[0 1.2*maxval])
                plot([0,0],[0, maxval],'k')
            end
    
    %                 plot([0,0],[0.5, 2],'k')
    %                 xlim(ax_tc,[0.5 3.5])
    %                 ylim([0.75, 1.5])
    
            xlim(ax_tc,[-tc_pre_sec, tc_post_sec])
            
    %                 ylabel('Norm FR')
            xlabel('Seconds since capture')                
            ylabel(ax_tc,'Firing Rate (Hz)')
            title_str = ['Anim: ',anim_fields{anim_i},'  Hunt Sesh: ',num2str(beh_i)];
            title(ax_tc,title_str);    
        end

    end

    drawnow    
    
    if Run_BOOT_ON == 1
        full_boot_percs{anim_i} = anim_boot_percs;
        full_boot_sigs{anim_i} = anim_boot_sigs;
    end

    full_cap_timecourses{anim_i} = anim_cap_timecourses;
    full_behFRs{anim_i} = anim_behFrs;
end


%% create bootstrap cell sig modulation charts
cell_sigs = cell(1,5);
for anim_i = 1:length(full_boot_sigs)
    anim_sigs = full_boot_sigs{anim_i};

    for beh_i = 1:length(anim_sigs)
        cell_sigs{beh_i} = [cell_sigs{beh_i}; anim_sigs{beh_i}];
    end
end

figure('Position',[36 301 836 307]); hold on;
full_pur_sigs = [];
full_cons_sigs = [];
for beh_i = 1:length(cell_sigs)
    ax_pie = subplot(2,length(cell_sigs),beh_i);
    
    pur_beh_sigs = cell_sigs{beh_i}(:,1);
    nonnan_inds = ~isnan(pur_beh_sigs);
    pur_beh_sigs = pur_beh_sigs(nonnan_inds);
    cons_beh_sigs = cell_sigs{beh_i}(nonnan_inds,2);
    pur_data_toplot = [sum(pur_beh_sigs==0), sum(pur_beh_sigs==1), sum(pur_beh_sigs==-1)];
    
    p1 = pie(pur_data_toplot);
    p1(1).FaceColor = [0.6 0.6 0.6];
    p1(3).FaceColor = c_pur*1.1;
    p1(5).FaceColor = c_pur*0.75;

    ax_pie2 = subplot(2,length(cell_sigs),beh_i+length(cell_sigs));
    cons_data_toplot = [sum(cons_beh_sigs==0), sum(cons_beh_sigs==1), sum(cons_beh_sigs==-1)];
    p2 = pie(cons_data_toplot);
    p2(1).FaceColor = [0.6 0.6 0.6];
    p2(3).FaceColor = c_cons;
    p2(5).FaceColor = c_cons*0.75;

    %     xlabel('Seconds since capture')                
%     ylabel(ax_tc,'Firing Rate (Norm to baseline)')
    title_str = ['Hunt Sesh: ',num2str(beh_i)];
    title(ax_pie,title_str);    
    
    full_pur_sigs = padmat(full_pur_sigs,cell_sigs{beh_i}(:,1),2);
    full_cons_sigs = padmat(full_cons_sigs,cell_sigs{beh_i}(:,2),2);
end

sgtitle(['All animal bootstrap sig (pursuit and cons.) ',num2str(num_boots),'iters ',' cell=',celltype])

% t = table(full_pur_sigs(:,1),full_pur_sigs(:,2),full_pur_sigs(:,3),...
%             'VariableNames',{'HS1','HS2','HS3'});
% Sesh = [1, 2, 3];
% % Fit the repeated measures model
% rmModel = fitrm(t, 'HS1-HS3~1', 'WithinDesign', Sesh);
% rmANOVA = ranova(rmModel) % general test for all conditions
% 
% % tests specific pairs
% [c] = multcompare(rmModel,'Time') % defaults to tukey-kramer critical value

% p_pur12 = c.pValue(1);
% p_pur13 = c.pValue(2);

% % repeat for cons
% t = table(full_cons_sigs(:,1),full_cons_sigs(:,2),full_cons_sigs(:,3),...
%             'VariableNames',{'HS1','HS2','HS3'});
% Sesh = [1, 2, 3];
% % Fit the repeated measures model
% rmModel = fitrm(t, 'HS1-HS3~1', 'WithinDesign', Sesh);
% rmANOVA = ranova(rmModel) % general test for all conditions
% 
% % tests specific pairs
% [c] = multcompare(rmModel,'Time') % defaults to tukey-kramer critical value
% 
% p_cons12 = c.pValue(1);
% p_cons13 = c.pValue(2);

full_pur_sigs_3col = full_pur_sigs(:,1:3);
full_pur_sigs_3col(sum(isnan(full_pur_sigs_3col),2)>0,:) = [];

pur_beh_sigs = full_pur_sigs_3col(:,1);
pur_data_toplot = [sum(pur_beh_sigs==0), sum(pur_beh_sigs==1), sum(pur_beh_sigs==-1)];

[p,tbl,stats] = friedman(full_pur_sigs_3col);
[c] = multcompare(stats);
p_pur12 = c(1,end);
p_pur13 = c(2,end);

full_cons_sigs_3col = full_cons_sigs(:,1:3);
full_cons_sigs_3col(sum(isnan(full_cons_sigs_3col),2)>0,:) = [];

[p,tbl,stats] = friedman(full_cons_sigs_3col);
[c] = multcompare(stats);
p_cons12 = c(1,end);
p_cons13 = c(2,end);

disp(['p-values  --  pur12=',num2str(p_pur12),get_asterisks_from_pval(p_pur12),...
    '  pur13=',num2str(p_pur13),get_asterisks_from_pval(p_pur13)])
disp(['p-values  --  cons12=',num2str(p_cons12),get_asterisks_from_pval(p_cons12),...
    '  cons13=',num2str(p_cons13),get_asterisks_from_pval(p_cons13)])
% 
% [pur_p12,h,stats] = signrank(full_pur_sigs(:,1),full_pur_sigs(:,2))
% [pur_p13,h,stats] = signrank(full_pur_sigs(:,1),full_pur_sigs(:,3))
% [pur_p23,h,stats] = signrank(full_pur_sigs(:,2),full_pur_sigs(:,3))
% 
% [cons_p12,h,stats] = signrank(full_cons_sigs(:,1),full_cons_sigs(:,2))
% [cons_p13,h,stats] = signrank(full_cons_sigs(:,1),full_cons_sigs(:,3))
% [cons_p23,h,stats] = signrank(full_cons_sigs(:,2),full_cons_sigs(:,3))


figure('Position',[14 588 875 257]); hold on;
select_inds = [1,2,3,5];
other_inds = [4,6,7,8];

full_pair_sums = [];
for beh_i = 1:length(cell_sigs)
    ax_pie = subplot(1,length(cell_sigs),beh_i);
    
    pur_beh_sigs = cell_sigs{beh_i}(:,1);
    cons_beh_sigs = cell_sigs{beh_i}(:,2);

    overlap_mat = [0,0,0,0,0,0,0,0]; 
    not_active_sum = 0;
    for i = 1:size(pur_beh_sigs,1)
        pair = [pur_beh_sigs(i), cons_beh_sigs(i)];

        if all(pair == [1,1])
            overlap_mat(1) = overlap_mat(1) + 1;
        elseif all(pair == [1,0])
            overlap_mat(2) = overlap_mat(2) + 1;
        elseif all(pair == [0,1])
            overlap_mat(3) = overlap_mat(3) + 1;
        elseif all(pair == [-1,0])
            overlap_mat(4) = overlap_mat(4) + 1;
        elseif all(pair == [0,-1])
            overlap_mat(5) = overlap_mat(5) + 1;
        elseif all(pair == [1,-1])
            overlap_mat(6) = overlap_mat(6) + 1;
        elseif all(pair == [-1,1])
            overlap_mat(7) = overlap_mat(7) + 1;
        elseif all(pair == [-1,-1])
            overlap_mat(8) = overlap_mat(8) + 1;
        else
            not_active_sum = not_active_sum + 1;
        end
    end
    

    select_pairs = overlap_mat(select_inds);
%     select_pairs(5) = sum(overlap_mat(other_inds)) + not_active_sum;
    select_pairs(5) = sum(overlap_mat(other_inds));

    full_pair_sums(beh_i,:) = select_pairs;

    p1 = pie(select_pairs);
    p1(1).FaceColor = c_both;
    p1(3).FaceColor = c_pur*1.1;
    p1(5).FaceColor = c_cons;
    p1(7).FaceColor = c_cons*0.75;
    p1(9).FaceColor = [0.6 0.6 0.6];

    title_str = ['Hunt Sesh: ',num2str(beh_i)];
    title(ax_pie,title_str);    

end

sgtitle(['All animal bootstrap sig (overlap) ',num2str(num_boots),'iters ',' cell=',celltype])

both_up_cols = (full_pur_sigs_3col == 1) & (full_cons_sigs_3col == 1);
pur_up_cols = (full_pur_sigs_3col == 1) & (full_cons_sigs_3col == 0);
cons_up_cols = (full_pur_sigs_3col == 0) & (full_cons_sigs_3col == 1);
cons_down_cols = (full_pur_sigs_3col == 0) & (full_cons_sigs_3col == -1);

[h,both_p,stats] = cochranqtest(both_up_cols);
[h,pur_p,stats] = cochranqtest(pur_up_cols);
[h,cons_up_p,stats] = cochranqtest(cons_up_cols);
[h,cons_down_p,stats] = cochranqtest(cons_down_cols);

disp(['p-values  --  both-up=',num2str(both_p),get_asterisks_from_pval(both_p)])
disp(['p-values  --  pur-up=',num2str(pur_p),get_asterisks_from_pval(pur_p)])
disp(['p-values  --  cons-up=',num2str(cons_up_p),get_asterisks_from_pval(cons_up_p)])
disp(['p-values  --  cons-down=',num2str(cons_down_p),get_asterisks_from_pval(cons_down_p)])





% for i = 1:length(cell_sigs{1})
%     pFS(i).HS1_sig = cell_sigs{1}(i,:);
%     pFS(i).HS2_sig = cell_sigs{2}(i,:);
%     pFS(i).HS3_sig = cell_sigs{3}(i,:);
% end

% create bootstrap cell sig modulation charts
cell_sigs = cell(1,5);
for anim_i = 1:length(full_boot_percs)
    anim_sigs = full_boot_percs{anim_i};

    for beh_i = 1:length(anim_sigs)
        cell_sigs{beh_i} = [cell_sigs{beh_i}; anim_sigs{beh_i}];
    end
end


fig_percs = figure('Position',[120 263 503 411]);
ax = gca;
x_pts = [ones(1,length(cell_sigs{1})); 2*ones(1,length(cell_sigs{1})); 3*ones(1,length(cell_sigs{1}))];
y_pts = [cell_sigs{1}(:,1)'; cell_sigs{2}(:,1)'; cell_sigs{3}(:,1)'];
plot(x_pts,y_pts,'Color',[0 0 0 0.15]);
hold('on')

plot([1 2 3],nanmean(y_pts,2),'k','LineWidth',2.5)

plot(ones(1,length(cell_sigs{1})),cell_sigs{1}(:,1),'go');
plot(2*ones(1,length(cell_sigs{1})),cell_sigs{2}(:,1),'go');
plot(3*ones(1,length(cell_sigs{1})),cell_sigs{3}(:,1),'go');            

xticks(1:3)
set(gca,'XTickLabel',{'Hunt1','Hunt2','Hunt3'},'XTickLabelRotation',45,'FontSize',10);
ylabel('Bootstrap percentile')
box(gca,'off')
xlim([0.5 3.5])
ylim([0,105])
set(gca,'YScale','log')

% maxval = max([ICI_FRs Pur_FRs Cons_FRs]);
% if maxval ~= 0
% ylim(ax,[0 1.2*maxval])
% end


%% Timecourses

sig_plot_only = 1; %0 = plot all cells, 1 = plot only sign. modulated
sig_plot_YES = 0;

% create cell sig structure for TCs
cell_sigs = cell(1,5);
for anim_i = 1:length(full_boot_sigs)
    anim_sigs = full_boot_sigs{anim_i};

    for beh_i = 1:length(anim_sigs)
        cell_sigs{beh_i} = [cell_sigs{beh_i}; anim_sigs{beh_i}];
    end
end

% create animal cell average timecourse
cell_timecourses = cell(1,5);
for anim_i = 1:length(full_cap_timecourses)
    anim_tc = full_cap_timecourses{anim_i};

    for beh_i = 1:length(anim_tc)
        cell_timecourses{beh_i} = [cell_timecourses{beh_i}; anim_tc{beh_i}];
    end

end

figure('Position',[18 265 1270 333]); hold on;
for beh_i = 1:length(cell_timecourses)
    ax_tc = subplot(1,length(cell_timecourses),beh_i);
    hold(ax_tc,'on')    
    
    all_cell_timecourses = cell_timecourses{beh_i};
    if sig_plot_only == 1
        pur_beh_sigs = cell_sigs{beh_i}(:,1);
        if sig_plot_YES == 1
            cell_tc_toplot = all_cell_timecourses(pur_beh_sigs==1,:);
        else
            cell_tc_toplot = all_cell_timecourses(pur_beh_sigs~=1,:);
        end
    else
        cell_tc_toplot = all_cell_timecourses;
    end
    
    plot([0,0],[0, 3],'color',[0 0 0 0.35])
    alpha = 0.3;
%     sh_color = [0.35 0.8 0.35];
    if sig_plot_YES == 1
        sh_color = c_pur;
    else
        sh_color = [0.6 0.6 0.6];
    end
    s1 = stdshade(cell_tc_toplot,alpha,sh_color,tc_time,1,'sem');    
    s1.LineWidth = 1.75 ;
%     maxval = max(nanmean(cell_timecourses{beh_i}));
%     if maxval ~= 0
%         ylim(ax_tc,[0 1.2*maxval])
%     end
    ylim([0.5 2.5])

    xlim(ax_tc,[-tc_pre_sec, tc_post_sec])
    r1 = refline(0,1);
    r1.LineWidth = 1.5;
    r1.LineStyle = '--';
    r1.Color = 'k';
    
%                 ylabel('Norm FR')
    if center_beh == 1
        xlabel('Sec from pursuit start')        
    elseif center_beh == 2
        xlabel('Sec from capture')        
    end

    ylabel(ax_tc,'Firing Rate (Norm to baseline)')
    title_str = ['Hunt Sesh: ',num2str(beh_i)];
    title(ax_tc,title_str);    

end

if center_beh == 1
    sgtitle(['FR timecourse around pursuit start ',num2str(tc_pre_sec) ,' pre sec ','cell=',...
        celltype, ' sigmod=',num2str(sig_plot_YES)])
elseif center_beh == 2
    sgtitle(['FR timecourse around capture ',num2str(tc_pre_sec) ,' pre sec ','cell=',...
        celltype, ' sigmod=',num2str(sig_plot_YES)])
end



%% Timecourse overlay

sig_plot_only = 1; %0 = plot all cells, 1 = plot only sign. modulated

% create cell sig structure for TCs
cell_sigs = cell(1,5);
for anim_i = 1:length(full_boot_sigs)
    anim_sigs = full_boot_sigs{anim_i};

    for beh_i = 1:length(anim_sigs)
        cell_sigs{beh_i} = [cell_sigs{beh_i}; anim_sigs{beh_i}];
    end
end

% create animal cell average timecourse
cell_timecourses = cell(1,5);
for anim_i = 1:length(full_cap_timecourses)
    anim_tc = full_cap_timecourses{anim_i};

    for beh_i = 1:length(anim_tc)
        cell_timecourses{beh_i} = [cell_timecourses{beh_i}; anim_tc{beh_i}];
    end

end

figure('Position',[18 265 1270 333]); hold on;
for beh_i = 1:length(cell_timecourses)
    ax_tc = subplot(1,length(cell_timecourses),beh_i);
    hold(ax_tc,'on')    
    
    all_cell_timecourses = cell_timecourses{beh_i};
    if sig_plot_only == 1
        pur_beh_sigs = cell_sigs{beh_i}(:,1);
        cell_tc_toplot_sig = all_cell_timecourses(pur_beh_sigs==1,:);
        cell_tc_toplot_notsig = all_cell_timecourses(pur_beh_sigs~=1,:);
    else
        cell_tc_toplot = all_cell_timecourses;
    end
    
    plot([0,0],[0, 3],'color',[0 0 0 0.35])
    alpha = 0.3;
%     sh_color = [0.35 0.8 0.35];
    sh_color = c_pur;
    sh_color = [0.6 0.6 0.6];
    s1 = stdshade(cell_tc_toplot_notsig,alpha,sh_color,tc_time,1,'sem');    
    s1.LineWidth = 1.75 ;

    sh_color = c_pur;
    s2 = stdshade(cell_tc_toplot_sig,alpha,sh_color,tc_time,1,'sem');    
    s2.LineWidth = 1.75 ;    
%     maxval = max(nanmean(cell_timecourses{beh_i}));
%     if maxval ~= 0
%         ylim(ax_tc,[0 1.2*maxval])
%     end
    ylim([0.5 2.5])

    xlim(ax_tc,[-tc_pre_sec, tc_post_sec])
    r1 = refline(0,1);
    r1.LineWidth = 1.5;
    r1.LineStyle = '--';
    r1.Color = 'k';
    
%                 ylabel('Norm FR')
    if center_beh == 1
        xlabel('Sec from pursuit start')        
    elseif center_beh == 2
        xlabel('Sec from capture')        
    end

%     legend([s1,s2],{'non-resp','resp'},'Location','northwest','FontSize',6)

    ylabel(ax_tc,'Firing Rate (Norm to baseline)')
    title_str = ['Hunt Sesh: ',num2str(beh_i)];
    title(ax_tc,title_str);    

end

if center_beh == 1
    sgtitle(['FR timecourse overlayed around pursuit start ',num2str(tc_pre_sec) ,' pre sec ','cell=',...
        celltype])
elseif center_beh == 2
    sgtitle(['FR timecourse overlayed around capture ',num2str(tc_pre_sec) ,' pre sec ','cell=',...
        celltype])
end





%% Pur Cons beh FR graphs

% create cell sig structure for TCs
cell_behFrs = cell(1,5);
for anim_i = 1:length(full_behFRs)
    anim_behFRs = full_behFRs{anim_i};

    for beh_i = 1:length(anim_behFRs)
        cell_behFrs{beh_i} = [cell_behFrs{beh_i}; anim_behFRs{beh_i}];
    end
end


figure('Position',[34 537 387 386]); hold on;
beh_ind = 1; %1 for pur, 2 for cons;
pur_colors = [c_pur; c_pur; c_pur];
alphas = [0.48,0.6,0.8];
full_data = [];
for beh_num = 1:3
    
    beh_data = cell_behFrs{beh_num}(:,beh_ind);
    beh_data(beh_data==Inf) = NaN;
    num_cells = sum(~isnan(beh_data));
    beh_mean = nanmean(beh_data);
    beh_err = nanstd(beh_data)/(sqrt(num_cells));
    beh_c = pur_colors(beh_num,:);

    beh_bar = bar(beh_num,beh_mean,0.6,'edgecolor',beh_c,'facecolor',beh_c,...
        'facealpha',alphas(beh_num),'linewidth',2);
    beh_err = errorbar(beh_num,beh_mean,beh_err,'linestyle','none','CapSize',15,...
        'color',beh_c,'linewidth',2);

    full_data = padmat(full_data,beh_data,2);
end

xlim([0.5 3.5])
xticks(1:3)
xticklabels({'HS1','HS2','HS3'})
title('Pursuit FR change')

if fr_change_calc == 2
    ylim([-10 65])
    ylabel('FR change Pur v ICI (%)')

end

t = table(full_data(:,1),full_data(:,2),full_data(:,3),...
            'VariableNames',{'HS1','HS2','HS3'});
Sesh = [1, 2, 3];
% Fit the repeated measures model
rmModel = fitrm(t, 'HS1-HS3~1', 'WithinDesign', Sesh);
rmANOVA = ranova(rmModel) % general test for all conditions

% tests specific pairs
[c] = multcompare(rmModel,'Time') % defaults to tukey-kramer critical value

p_pur12 = c.pValue(1);
p_pur13 = c.pValue(2);
p_pur23 = c.pValue(4);



% [h,p,ci,stats] = ttest2(full_data(:,1),full_data(:,2))
% [h,p,ci,stats] = ttest2(full_data(:,1),full_data(:,3))
% [h,p,ci,stats] = ttest2(full_data(:,2),full_data(:,3))


figure('Position',[34 537 387 386]); hold on;
beh_ind = 2; %1 for pur, 2 for cons;
pur_colors = [c_cons; c_cons; c_cons].*0.95;
alphas = [0.5,0.67,0.81];
full_data = [];
for beh_num = 1:3
    
    beh_data = cell_behFrs{beh_num}(:,beh_ind);
    beh_data(beh_data==Inf) = NaN;
    num_cells = sum(~isnan(beh_data));
    beh_mean = nanmean(beh_data);
    beh_err = nanstd(beh_data)/(sqrt(num_cells));
    beh_c = pur_colors(beh_num,:);

    beh_bar = bar(beh_num,beh_mean,0.6,'edgecolor',beh_c,'facecolor',beh_c,...
        'facealpha',alphas(beh_num),'linewidth',2);
    beh_err = errorbar(beh_num,beh_mean,beh_err,'linestyle','none','CapSize',15,...
        'color',beh_c,'linewidth',2);

    full_data = padmat(full_data,beh_data,2);
end

xlim([0.5 3.5])
xticks(1:3)
xticklabels({'HS1','HS2','HS3'})
title('Consumption FR change')

if fr_change_calc == 2
    ylim([-10 10])
    ylim([-15 10])
    ylabel('FR change Consumption v ICI (%)')

end

t = table(full_data(:,1),full_data(:,2),full_data(:,3),...
            'VariableNames',{'HS1','HS2','HS3'});
Sesh = [1, 2, 3];
% Fit the repeated measures model
rmModel = fitrm(t, 'HS1-HS3~1', 'WithinDesign', Sesh);
rmANOVA = ranova(rmModel) % general test for all conditions

% tests specific pairs
[c] = multcompare(rmModel,'Time') % defaults to tukey-kramer critical value

p_cons12 = c.pValue(1);
p_cons13 = c.pValue(2);
p_cons23 = c.pValue(4);


disp(['p-values  --  pur12=',num2str(p_pur12),get_asterisks_from_pval(p_pur12),...
    '  pur13=',num2str(p_pur13),get_asterisks_from_pval(p_pur13),...
    '  pur23=',num2str(p_pur23),get_asterisks_from_pval(p_pur23)])

disp(['p-values  --  cons12=',num2str(p_cons12),get_asterisks_from_pval(p_cons12),...
    '  cons13=',num2str(p_cons13),get_asterisks_from_pval(p_cons13),...
    '  cons23=',num2str(p_cons23),get_asterisks_from_pval(p_cons23)])


% [h,p,ci,stats] = ttest2(full_data(:,1),full_data(:,2))
% [h,p,ci,stats] = ttest2(full_data(:,1),full_data(:,3))
% [h,p,ci,stats] = ttest2(full_data(:,2),full_data(:,3))

