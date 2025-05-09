%% DREADDs PIPELINE for Leman et al., 2025
% Pipeline for analyzing behavioral performance of Saline vs. XPR01595
% treated animals as part of Figure 7 in Leman et al., 2025


xp_code_dir = pwd; % run pipeline in DREADDs directory containing code and the folder "animal_data"
% uigetdir(); % Choose directory containing DREADDs code
addpath(genpath(xp_code_dir))

%% XPRO

xpro_times = dir('animal_data/xpro/*.csv');
xpro = struct2cell(xpro_times);

%DREADDs Hunt Calc
[xpro_t2c_mat,  xpro_lat_mat,  xpro_ad_mat] = DLC_huntcalc(xpro, xpro_times);

%% CTL

ctl_times = dir('animal_data/ctl/*.csv');
ctl = struct2cell(ctl_times);

[ctl_t2c_mat,ctl_lat_mat, ctl_ad_mat] = DLC_huntcalc(ctl, ctl_times);

%% Compare Conditions and Generate Figures

cd(xp_code_dir)
type = 1; % use 1 for day averages
exp = 1; % 1 = time to capture, 2 = latency to attack, 3 = pursuit duration
hs = 4; % use 0 if for no plotting, or 1-3 for those sessions alone, or 4 for all
sess = 5; % use 5 for xpro (single day learning (3) + probe day (2))
if exp == 1
    ctl_mat = ctl_t2c_mat;
    xpro_mat = xpro_t2c_mat;

elseif exp == 2 
    ctl_mat = ctl_lat_mat;
    xpro_mat = xpro_lat_mat;
elseif exp == 3
    ctl_mat = ctl_ad_mat;
    xpro_mat = xpro_ad_mat;
end 

DLC_huntplot_xp(ctl_mat(:,1:sess), xpro_mat(:,1:sess), type, exp, hs)


