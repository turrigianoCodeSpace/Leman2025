% plot both wake and sleep delta FR post hunt
%
% Written by Brian Cary
%

% Run both extwake and extsleep analysis to created stored data workspaces

mfilePath = mfilename('fullpath');
if isempty(mfilePath)
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
cd(fileparts(mfilePath))

% Path to processed data structure including spike times, features, etc.
hunt_path = "Z:\DPL\PROJECTS\StateCoding\Hunt_MasterStrct_StableBase.mat";
sp_filepath = split(hunt_path,'\');
sp_path2 = split(sp_filepath(end),'.mat');
savedata_dir = [fileparts(mfilePath), '\stored_data'];
savename = [char(sp_path2{1}),'_extwake_WS.mat'];
hunt_wake_path = [savedata_dir,filesep,savename];

savename = [char(sp_path2{1}),'_extsleep_WS.mat'];
hunt_sleep_path = [savedata_dir,filesep,savename];

cont_path = "Z:\DPL\PROJECTS\StateCoding\Cont_MasterStrct_StableBase.mat";
sp_filepath = split(cont_path,'\');
sp_path2 = split(sp_filepath(end),'.mat');
savedata_dir = [fileparts(mfilePath), '\stored_data'];
savename = [char(sp_path2{1}),'_extwake_WS.mat'];
cont_wake_path = [savedata_dir,filesep,savename];

savename = [char(sp_path2{1}),'_extsleep_WS.mat'];
cont_sleep_path = [savedata_dir,filesep,savename];

% RSU meant6
cont_wake_ws = load(cont_wake_path);
cont_sleep_ws = load(cont_sleep_path);
hunt_wake_ws = load(hunt_wake_path);
hunt_sleep_ws = load(hunt_sleep_path);


rsu = 1;
mean_t = 6;
test_type = ' '; % ttest or signrank

% colors for plotting
c_rem   = [25 181 149]./255;
c_nrem  = [131 49 146]./255;
c_aw    = [201 28 101]./255;
c_qw    = [247 148 41]./255;

c_rem_cont   = c_rem*0.75;
c_nrem_cont  = c_nrem*0.75;
c_aw_cont    = c_aw*0.75;
c_qw_cont    = c_qw*0.75;

set(0,'defaultFigureUnit','pixels');
set(0,'defaultFigurePosition',[82 584 358 361]);

%% Fig  - plotting the first-last change
% SETTINGS FOR Z-SCORE
%% HUNT

mylims = 1;

AW_delta_nonan = hunt_wake_ws.AW_delta_nonan;
QW_delta_nonan = hunt_wake_ws.QW_delta_nonan;
NREM_delta_nonan = hunt_sleep_ws.NREM_delta_nonan;

AW_mean = nanmean(AW_delta_nonan);
AW_sem = std(AW_delta_nonan,0,'omitnan') / sqrt(sum(~isnan(AW_delta_nonan)) - 1);

QW_mean = nanmean(QW_delta_nonan);
QW_sem = std(QW_delta_nonan,0,'omitnan') / sqrt(sum(~isnan(QW_delta_nonan)) - 1);

NREM_mean = nanmean(NREM_delta_nonan);
NREM_sem = std(NREM_delta_nonan,0,'omitnan') / sqrt(sum(~isnan(NREM_delta_nonan)) - 1);

% REM_mean = nanmean(REM_delta_nonan);
% REM_sem = std(REM_delta_nonan,0,'omitnan') / sqrt(sum(~isnan(REM_delta_nonan)) - 1);

if strcmp(test_type,'ttest')
    [~,p_aw] = ttest(AW_delta_nonan);
    [~,p_qw] = ttest(QW_delta_nonan);
    [~,p_nrem] = ttest(NREM_delta_nonan);
% 
else
    [p_aw,h,stats] = signrank(AW_delta_nonan);
    [p_qw,h,stats] = signrank(QW_delta_nonan);
    [p_nrem,h,stats] = signrank(NREM_delta_nonan);
end

p_aw = p_aw*3;
p_qw = p_qw*3;
p_nrem = p_nrem*3;
[a_aw,fsz_aw] = get_asterisks_from_pval(p_aw);
[a_qw,fsz_qw] = get_asterisks_from_pval(p_qw);
[a_nrem,fsz_nrem] = get_asterisks_from_pval(p_nrem);

disp(['p-values  --  AW=',num2str(p_aw),'  QW=',num2str(p_qw),'  NREM=',num2str(p_nrem)])

dfig = figure();
set(dfig,'position',[102 528 348 399],'units','pixels');
box off
hold on;
bw = .5;
csz = 20;
mfc1 = c_aw;
mfc2 = c_qw;
ec1 = 'none';
ec2 = 'none';

aw_bar = bar(1,AW_mean,bw,'edgecolor',ec1,'facecolor',mfc1,'linewidth',2);
aw_err = errorbar(1,AW_mean,AW_sem,'linestyle','none','capsize',csz,...
    'color',c_aw,'linewidth',2);
qw_bar = bar(2,QW_mean,bw,'edgecolor',ec2,'facecolor',mfc2,'linewidth',2);
qw_err = errorbar(2,QW_mean,QW_sem,'linestyle','none','capsize',csz,...
    'color',c_qw,'linewidth',2);
nrem_bar = bar(3,NREM_mean,bw,'edgecolor',ec2,'facecolor',c_nrem,'linewidth',2);
nrem_err = errorbar(3,NREM_mean,NREM_sem,'linestyle','none','capsize',csz,...
    'color',c_nrem,'linewidth',2);

if mylims
    if mean_t == 4
        xl = [.5 3.5];
        yl = [-0.1 0.15];
        yt = -0.1:.05:0.15;
    elseif mean_t == 6
        if rsu == 1
            xl = [.5 3.5];
            yl = [-10 35];
            yt = [-10:10:35];
        else
            xl = [.5 3.5];
            yl = [-10 25];
            yt = [-10:10:25];
        end
    end

else
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    yt = get(gca,'ytick');
end

if mean_t == 4
    ystar_ht = 0.12;
elseif mean_t == 6
    if rsu == 1
        ystar_ht = 28;
    else
        ystar_ht = 18;
    end
end

text(1,ystar_ht,a_aw,'fontsize',18,'Rotation',45);
text(2,ystar_ht,a_qw,'fontsize',18,'Rotation',45);
text(3,ystar_ht,a_nrem,'fontsize',18,'rotation',45); 
set(gca,'xlim',xl,'ylim',yl,'xtick',[1 2 3],'xticklabel',{'Active','Quiet','NREM'},...
    'ytick',yt);
if mean_t == 4
    ylabel('Firing rate (z)','fontsize',18);
elseif mean_t == 6
    ylabel({'Firing rate change','(% change from first QW epoch)'},'fontsize',14);
end
title('Hunt - Ext. Sleep/Wake delta FR by state','fontsize',12);


x = (AW_delta_nonan-nanmean(AW_delta_nonan))/nanstd(AW_delta_nonan);
%         log_dFR = log(AW_dFR_nonan);
%         x = (log_dFR-nanmean(log_dFR))/nanstd(log_dFR);
[h,p] = kstest(x);
figure;hist(x,100);title('z-score dFR AW histogram')
disp(['AW dFR normality rejection p-val: ',num2str(p)])

%% CONTROL
AW_delta_nonan = cont_wake_ws.AW_delta_nonan;
QW_delta_nonan = cont_wake_ws.QW_delta_nonan;
NREM_delta_nonan = cont_sleep_ws.NREM_delta_nonan;

AW_mean = nanmean(AW_delta_nonan);
AW_sem = std(AW_delta_nonan,0,'omitnan') / sqrt(sum(~isnan(AW_delta_nonan)) - 1);

QW_mean = nanmean(QW_delta_nonan);
QW_sem = std(QW_delta_nonan,0,'omitnan') / sqrt(sum(~isnan(QW_delta_nonan)) - 1);

NREM_mean = nanmean(NREM_delta_nonan);
NREM_sem = std(NREM_delta_nonan,0,'omitnan') / sqrt(sum(~isnan(NREM_delta_nonan)) - 1);

% REM_mean = nanmean(REM_delta_nonan);
% REM_sem = std(REM_delta_nonan,0,'omitnan') / sqrt(sum(~isnan(REM_delta_nonan)) - 1);


if strcmp(test_type,'ttest')
    [~,p_aw] = ttest(AW_delta_nonan);
    [~,p_qw] = ttest(QW_delta_nonan);
    [~,p_nrem] = ttest(NREM_delta_nonan);
% 
else
    [p_aw,h,stats] = signrank(AW_delta_nonan);
    [p_qw,h,stats] = signrank(QW_delta_nonan);
    [p_nrem,h,stats] = signrank(NREM_delta_nonan);
end

p_aw = p_aw*3;
p_qw = p_qw*3;
p_nrem = p_nrem*3;
[a_aw,fsz_aw] = get_asterisks_from_pval(p_aw);
[a_qw,fsz_qw] = get_asterisks_from_pval(p_qw);
[a_nrem,fsz_nrem] = get_asterisks_from_pval(p_nrem);

disp(['p-values  --  AW=',num2str(p_aw),'  QW=',num2str(p_qw),'  NREM=',num2str(p_nrem)])

dfig = figure();
set(dfig,'position',[102 528 348 399],'units','pixels');
box off
hold on;
bw = .5;
csz = 20;
mfc1 = c_aw_cont;
mfc2 = c_qw_cont;
ec1 = 'none';
ec2 = 'none';

aw_bar = bar(1,AW_mean,bw,'edgecolor',ec1,'facecolor',mfc1,'linewidth',2);
aw_err = errorbar(1,AW_mean,AW_sem,'linestyle','none','capsize',csz,...
    'color',c_aw_cont,'linewidth',2);
qw_bar = bar(2,QW_mean,bw,'edgecolor',ec2,'facecolor',mfc2,'linewidth',2);
qw_err = errorbar(2,QW_mean,QW_sem,'linestyle','none','capsize',csz,...
    'color',c_qw_cont,'linewidth',2);
nrem_bar = bar(3,NREM_mean,bw,'edgecolor',ec2,'facecolor',c_nrem_cont,'linewidth',2);
nrem_err = errorbar(3,NREM_mean,NREM_sem,'linestyle','none','capsize',csz,...
    'color',c_nrem_cont,'linewidth',2);

if mylims
    if mean_t == 4
        xl = [.5 3.5];
        yl = [-0.1 0.15];
        yt = -0.1:.05:0.15;
    elseif mean_t == 6
        if rsu == 1
            xl = [.5 3.5];
            yl = [-10 35];
            yt = [-10:10:35];
        else
            xl = [.5 3.5];
            yl = [-10 25];
            yt = [-10:10:25];
        end
    end

else
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    yt = get(gca,'ytick');
end

if mean_t == 4
    ystar_ht = 0.12;
elseif mean_t == 6
    if rsu == 1
        ystar_ht = 20;
    else
        ystar_ht = 18;
    end
end

text(1,ystar_ht,a_aw,'fontsize',18,'Rotation',45);
text(2,ystar_ht,a_qw,'fontsize',18,'Rotation',45);
text(3,ystar_ht,a_nrem,'fontsize',18,'rotation',45); 

set(gca,'xlim',xl,'ylim',yl,'xtick',[1 2 3],'xticklabel',{'Active','Quiet','NREM'},...
    'ytick',yt);
ylabel('Firing rate change (z)','fontsize',16);
title('CONTROL - Ext. Sleep/Wake delta FR by state','fontsize',12);



%% POST HUNT GROUPED DELTAS
%% HUNT Z-SCORE
% load in vars
AW_delta_nonan = hunt_wake_ws.AW_delta_nonan;
QW_delta_nonan = hunt_wake_ws.QW_delta_nonan;
NREM_delta_nonan = hunt_sleep_ws.NREM_delta_nonan;
REM_delta_nonan = hunt_sleep_ws.REM_delta_nonan;

wake_starts_nonan_AW = hunt_wake_ws.wake_starts_nonan_AW;
wake_starts_nonan_QW = hunt_wake_ws.wake_starts_nonan_QW;
sleep_starts_nonan_NREM = hunt_sleep_ws.sleep_starts_nonan_NREM;
sleep_starts_nonan_REM = hunt_sleep_ws.sleep_starts_nonan_REM;

wakedur_nonan_AW = hunt_wake_ws.wakedur_nonan_AW;
wakedur_nonan_QW = hunt_wake_ws.wakedur_nonan_QW;
sleepdur_nonan_NREM = hunt_sleep_ws.sleepdur_nonan_NREM;
sleepdur_nonan_REM = hunt_sleep_ws.sleepdur_nonan_REM;

post_hunt_hr = hunt_wake_ws.post_hunt_hr;

wake_starts_AW = wake_starts_nonan_AW./3600;
wake_starts_QW = wake_starts_nonan_QW./3600;

wake_ends_AW = wake_starts_AW + wakedur_nonan_AW/3600;
wake_ends_QW = wake_starts_QW + wakedur_nonan_QW/3600;

% calculate light vs dark circ data bins
AW_this_group = {};
QW_this_group = {};
AW_group_means = [];
AW_group_sem = [];
QW_group_means = [];
QW_group_sem = [];
group_means_offT = [];

group_hrs = 12;
group_split = post_hunt_hr : group_hrs : 6.5*24;
ngroups = numel(group_split)-1;
for uu = 1:ngroups
    
    g0 = group_split(uu);
    g1 = group_split(uu+1);
    AW_this_group_idx = find(wake_starts_AW >= g0 & wake_ends_AW < g1);
    QW_this_group_idx = find(wake_starts_QW >= g0 & wake_ends_QW < g1);
    AW_this_group{uu} = AW_delta_nonan(AW_this_group_idx);
    QW_this_group{uu}= QW_delta_nonan(QW_this_group_idx);
    AW_group_means(uu) = nanmean(AW_this_group{uu});
    AW_group_sem(uu) = std(AW_this_group{uu},0,'omitnan') / sqrt(numel(AW_this_group{uu})-1);
    QW_group_means(uu) = nanmean(QW_this_group{uu});
    QW_group_sem(uu) = std(QW_this_group{uu},0,'omitnan') / sqrt(numel(QW_this_group{uu})-1);
    group_means_offT(uu) = mean([g0,g1]);
end

sleep_starts_REM = sleep_starts_nonan_REM./3600;
sleep_starts_NREM = sleep_starts_nonan_NREM./3600;

sleep_ends_REM = sleep_starts_REM + sleepdur_nonan_REM/3600;
sleep_ends_NREM = sleep_starts_NREM + sleepdur_nonan_NREM/3600;

% calculate light vs dark circ data bins
REM_this_group = {};
NREM_this_group = {};
REM_group_means = [];
REM_group_sem = [];
NREM_group_means = [];
NREM_group_sem = [];
group_means_offT = [];

group_hrs = 12;
group_split = post_hunt_hr : group_hrs : 6.5*24;
ngroups = numel(group_split)-1;
for uu = 1:ngroups
    
    g0 = group_split(uu);
    g1 = group_split(uu+1);
    REM_this_group_idx = find(sleep_starts_REM >= g0 & sleep_ends_REM < g1);
    NREM_this_group_idx = find(sleep_starts_NREM >= g0 & sleep_ends_NREM < g1);
    REM_this_group{uu} = REM_delta_nonan(REM_this_group_idx);
    NREM_this_group{uu}= NREM_delta_nonan(NREM_this_group_idx);
    REM_group_means(uu) = nanmean(REM_this_group{uu});
    REM_group_sem(uu) = std(REM_this_group{uu},0,'omitnan') / sqrt(numel(REM_this_group{uu})-1);
    NREM_group_means(uu) = nanmean(NREM_this_group{uu});
    NREM_group_sem(uu) = std(NREM_this_group{uu},0,'omitnan') / sqrt(numel(NREM_this_group{uu})-1);
    group_means_offT(uu) = mean([g0,g1]);
end

%% PLOT
dfig = figure();
set(dfig,'position',[48 158 650 361]);
box off
hold on;
bw = .4;
csz = 12;
ec1 = 'none';
ec2 = 'none';    
mylims = 0;
for per_num = 1:length(group_split)-1

    if strcmp(test_type,'ttest')
        [~,p_aw] = ttest(AW_this_group{per_num});
        [~,p_qw] = ttest(QW_this_group{per_num});
        [~,p_nrem] = ttest(NREM_this_group{per_num});
    else
        [p_aw,h,stats] = signrank(AW_this_group{per_num});
        [p_qw,h,stats] = signrank(QW_this_group{per_num});
        [p_nrem,h,stats] = signrank(NREM_this_group{per_num});
    end    
    
    p_aw = p_aw*3;
    p_qw = p_qw*3;
    p_nrem = p_nrem*3;
    [a_aw,fsz_aw] = get_asterisks_from_pval(p_aw);
    [a_qw,fsz_qw] = get_asterisks_from_pval(p_qw);
    [a_nrem,fsz_nrem] = get_asterisks_from_pval(p_nrem);
    
    disp(['group=',num2str(per_num),'  p-values  --  AW=',num2str(p_aw),'  QW=',num2str(p_qw),'  NREM=',num2str(p_nrem)])
    
    bar_x1 = (per_num-1)*2 + 1;
    bar_x2 = (per_num-1)*2 + 1.5;
    bar_x3 = (per_num-1)*2 + 2;
    bar_x4 = (per_num-1)*2 + 2.5;

    aw_bar = bar(bar_x1,AW_group_means(per_num),bw,'edgecolor',ec1,'facecolor',c_aw,'linewidth',2);
    aw_err = errorbar(bar_x1,AW_group_means(per_num),AW_group_sem(per_num),'linestyle','none','capsize',csz,...
        'color',c_aw,'linewidth',2);
    qw_bar = bar(bar_x2,QW_group_means(per_num),bw,'edgecolor',ec2,'facecolor',c_qw,'linewidth',2);
    qw_err = errorbar(bar_x2,QW_group_means(per_num),QW_group_sem(per_num),'linestyle','none','capsize',csz,...
        'color',c_qw,'linewidth',2);
    nrem_bar = bar(bar_x3,NREM_group_means(per_num),bw,'edgecolor',ec2,'facecolor',c_nrem,'linewidth',2);
    nrem_err = errorbar(bar_x3,NREM_group_means(per_num),NREM_group_sem(per_num),'linestyle','none','capsize',csz,...
        'color',c_nrem,'linewidth',2);
%     REM_bar = bar(bar_x4,REM_group_means(per_num),bw,'edgecolor',ec1,'facecolor',c_rem,'linewidth',2);
%     REM_err = errorbar(bar_x4,REM_group_means(per_num),REM_group_sem(per_num),'linestyle','none','capsize',csz,...
%         'color',c_rem,'linewidth',2);

%     text(bar_x1,0.2,sprintf('p = %.4f',p_aw),'fontsize',10,'Rotation',45);
%     text(bar_x2,0.2,sprintf('p = %.4f',p_qw),'fontsize',10,'Rotation',45);
%     text(bar_x4,0.05,sprintf('p = %.4f',p_rem),'fontsize',10,'rotation',45);
%     text(bar_x3,0.05,sprintf('p = %.4f',p_nrem),'fontsize',10,'rotation',45);   
        
    if mean_t == 4
        ystar_ht = 0.16;
    elseif mean_t == 6
        if rsu == 1
            ystar_ht = 35;
        else
            ystar_ht = 25;
        end
    end

    text(bar_x1,ystar_ht,a_aw,'fontsize',18,'Rotation',45);
    text(bar_x2,ystar_ht,a_qw,'fontsize',18,'Rotation',45);
    text(bar_x3,ystar_ht,a_nrem,'fontsize',18,'rotation',45); 
%     text(bar_x4,ystar_ht,a_rem,'fontsize',18,'rotation',45);

%     ylabel('Firing rate change (z)','fontsize',16);

end

post_hunt_hrs = group_split+group_hrs/2 - post_hunt_hr;
set(gca,'xtick',1.5:2:length(group_split)*2,'xticklabel',post_hunt_hrs)
xlabel('Hours post hunt')
if mean_t == 4
    ylabel('Firing rate (z)','fontsize',18);
elseif mean_t == 6
    ylabel({'Firing rate change','(% change from first QW epoch)'},'fontsize',14);
end

title('Hunt - FR change by state post hunt')
if mean_t == 4
    ylim([-.1,0.2]);
elseif mean_t == 6
    ylim([-15, 35]);
    if rsu == 1
        ylim([-10, 55]);
    else
        ylim([-15, 35]);
    end    
end
xlim([0.5, 6.5])

% ylim([-.15,0.3])



%% Older versions

dfig = figure();
set(dfig,'position',[0.0510 0.4944 0.3969 0.3407]);
box off
hold on;
bw = .3;
csz = 12;
ec1 = 'none';
ec2 = 'none';    
mylims = 0;
for per_num = 1:length(group_split)-1
    %light group == 1
%         per_num = 1;
    
    if ~isempty(AW_this_group{per_num})
        [~,p_aw] = ttest(AW_this_group{per_num});
        [a_aw,fsz_aw] = get_asterisks_from_pval(p_aw);
    else
        p_aw = 'nan';
    end
    if ~isempty(QW_this_group{per_num})
        [~,p_qw] = ttest(QW_this_group{per_num});
        [a_qw,fsz_qw] = get_asterisks_from_pval(p_qw);
    else
        p_qw = 'nan';
    end
    if ~isempty(REM_this_group{per_num})
        [~,p_rem] = ttest(REM_this_group{per_num});
        [a_rem,fsz_rem] = get_asterisks_from_pval(p_aw);
    else
        p_rem = 'nan';
    end
    if ~isempty(NREM_this_group{per_num})
        [~,p_nrem] = ttest(NREM_this_group{per_num});
        [a_nrem,fsz_nrem] = get_asterisks_from_pval(p_nrem);
    else
        p_nrem = 'nan';
    end

    
    bar_x1 = (per_num-1)*2 + 1;
    bar_x2 = (per_num-1)*2 + 1.5;
    bar_x3 = (per_num-1)*2 + 2;
    bar_x4 = (per_num-1)*2 + 2.5;

    aw_bar = bar(bar_x1,AW_group_means(per_num),bw,'edgecolor',ec1,'facecolor',c_aw,'linewidth',2);
    aw_err = errorbar(bar_x1,AW_group_means(per_num),AW_group_sem(per_num),'linestyle','none','capsize',csz,...
        'color',c_aw,'linewidth',2);
    qw_bar = bar(bar_x2,QW_group_means(per_num),bw,'edgecolor',ec2,'facecolor',c_qw,'linewidth',2);
    qw_err = errorbar(bar_x2,QW_group_means(per_num),QW_group_sem(per_num),'linestyle','none','capsize',csz,...
        'color',c_qw,'linewidth',2);
    nrem_bar = bar(bar_x3,NREM_group_means(per_num),bw,'edgecolor',ec2,'facecolor',c_nrem,'linewidth',2);
    nrem_err = errorbar(bar_x3,NREM_group_means(per_num),NREM_group_sem(per_num),'linestyle','none','capsize',csz,...
        'color',c_nrem,'linewidth',2);
%     REM_bar = bar(bar_x4,REM_group_means(per_num),bw,'edgecolor',ec1,'facecolor',c_rem,'linewidth',2);
%     REM_err = errorbar(bar_x4,REM_group_means(per_num),REM_group_sem(per_num),'linestyle','none','capsize',csz,...
%         'color',c_rem,'linewidth',2);

%     text(bar_x1,0.2,sprintf('p = %.4f',p_aw),'fontsize',10,'Rotation',45);
%     text(bar_x2,0.2,sprintf('p = %.4f',p_qw),'fontsize',10,'Rotation',45);
%     text(bar_x4,0.05,sprintf('p = %.4f',p_rem),'fontsize',10,'rotation',45);
%     text(bar_x3,0.05,sprintf('p = %.4f',p_nrem),'fontsize',10,'rotation',45);   

    text(bar_x1,0.2,a_aw,'fontsize',10,'Rotation',45);
    text(bar_x2,0.2,a_qw,'fontsize',10,'Rotation',45);
    text(bar_x3,0.2,a_nrem,'fontsize',10,'rotation',45); 
%     text(bar_x4,0.2,a_rem,'fontsize',10,'rotation',45);

    ylabel('Firing rate change (z)','fontsize',16);

end

post_hunt_hrs = group_split+group_hrs/2 - post_hunt_hr;
set(gca,'xtick',1.5:2:length(group_split)*2,'xticklabel',post_hunt_hrs)
xlabel('Hours post hunt')
title('FR change by state post hunt')
ylim([-.1,0.3])


% meant = 5

dfig = figure();
set(dfig,'position',[0.0510 0.4944 0.3969 0.3407]);
box off
hold on;
bw = .3;
csz = 12;
ec1 = 'none';
ec2 = 'none';    
mylims = 0;
for per_num = 1:length(group_split)-1
    %light group == 1
%         per_num = 1;
    
    if ~isempty(AW_this_group{per_num})
        [~,p_aw] = ttest(AW_this_group{per_num}-1);
        [a_aw,fsz_aw] = get_asterisks_from_pval(p_aw);
    else
        p_aw = 'nan';
    end
    if ~isempty(QW_this_group{per_num})
        [~,p_qw] = ttest(QW_this_group{per_num}-1);
        [a_qw,fsz_qw] = get_asterisks_from_pval(p_qw);
    else
        p_qw = 'nan';
    end
    if ~isempty(REM_this_group{per_num})
        [~,p_rem] = ttest(REM_this_group{per_num}-1);
        [a_rem,fsz_rem] = get_asterisks_from_pval(p_aw);
    else
        p_rem = 'nan';
    end
    if ~isempty(NREM_this_group{per_num})
        [~,p_nrem] = ttest(NREM_this_group{per_num}-1);
        [a_nrem,fsz_nrem] = get_asterisks_from_pval(p_nrem);
    else
        p_nrem = 'nan';
    end

    
    bar_x1 = (per_num-1)*2 + 1;
    bar_x2 = (per_num-1)*2 + 1.5;
    bar_x3 = (per_num-1)*2 + 2;
    bar_x4 = (per_num-1)*2 + 2.5;

    aw_bar = bar(bar_x1,AW_group_means(per_num),bw,'edgecolor',ec1,'facecolor',c_aw,'linewidth',2);
    aw_err = errorbar(bar_x1,AW_group_means(per_num),AW_group_sem(per_num),'linestyle','none','capsize',csz,...
        'color',c_aw,'linewidth',2);
    qw_bar = bar(bar_x2,QW_group_means(per_num),bw,'edgecolor',ec2,'facecolor',c_qw,'linewidth',2);
    qw_err = errorbar(bar_x2,QW_group_means(per_num),QW_group_sem(per_num),'linestyle','none','capsize',csz,...
        'color',c_qw,'linewidth',2);
    nrem_bar = bar(bar_x3,NREM_group_means(per_num),bw,'edgecolor',ec2,'facecolor',c_nrem,'linewidth',2);
    nrem_err = errorbar(bar_x3,NREM_group_means(per_num),NREM_group_sem(per_num),'linestyle','none','capsize',csz,...
        'color',c_nrem,'linewidth',2);
%     REM_bar = bar(bar_x4,REM_group_means(per_num),bw,'edgecolor',ec1,'facecolor',c_rem,'linewidth',2);
%     REM_err = errorbar(bar_x4,REM_group_means(per_num),REM_group_sem(per_num),'linestyle','none','capsize',csz,...
%         'color',c_rem,'linewidth',2);

%     text(bar_x1,0.2,sprintf('p = %.4f',p_aw),'fontsize',10,'Rotation',45);
%     text(bar_x2,0.2,sprintf('p = %.4f',p_qw),'fontsize',10,'Rotation',45);
%     text(bar_x4,0.05,sprintf('p = %.4f',p_rem),'fontsize',10,'rotation',45);
%     text(bar_x3,0.05,sprintf('p = %.4f',p_nrem),'fontsize',10,'rotation',45);   

    text(bar_x1,1.35,a_aw,'fontsize',10,'Rotation',45);
    text(bar_x2,1.35,a_qw,'fontsize',10,'Rotation',45);
    text(bar_x3,1.35,a_nrem,'fontsize',10,'rotation',45); 
%     text(bar_x4,0.2,a_rem,'fontsize',10,'rotation',45);

    ylabel({['Firing rate change'],['(last / first epoch)']},'fontsize',16);

end

post_hunt_hrs = group_split+group_hrs/2 - post_hunt_hr;
set(gca,'xtick',1.5:2:length(group_split)*2,'xticklabel',post_hunt_hrs)
xlabel('Hours post hunt')
title('FR change by state post hunt')

ylim([0.8 1.6])
yticks(0.8:0.2:1.6)

ref_l = refline(0,1);
ref_l.LineStyle = '--';
ref_l.LineWidth = 1.5;
ref_l.Color = 'k';             
