%% DREADDs PIPELINE for Leman et al., 2025
% Pipeline for analyzing behavioral performance of PV-Cre+ and PV-Cre-
% negative animals as part of Figure 1 in Leman et al., 2025


dr_code_dir = pwd; % run pipeline in DREADDs directory containing code and the folder "animal_data"
% uigetdir(); % Choose directory containing DREADDs code
addpath(genpath(dr_code_dir))

%% Cre+

crep_times = dir('animal_data/cre+/*.csv');
crep = struct2cell(crep_times);

%DREADDs Hunt Calc
[dr_t2c_mat,  dr_lat_mat,  dr_ad_mat] = DLC_huntcalc(crep, crep_times);

%% Cre-

crem_times = dir('animal_data/cre-/*.csv');
crem = struct2cell(crem_times);

[ctl_t2c_mat,ctl_lat_mat, ctl_ad_mat] = DLC_huntcalc(crem, crem_times);

%% Compare Conditions and Generate Figures

cd(dr_code_dir)
type = 2; % use 1 for day averages, use 2 for session averages
exp = 1; % 1 = time to capture, 2 = latency to attack, 3 = pursuit duration
hs = 4; % use 0 if for no plotting, or 1-3 for those sessions alone, or 4 for all
sess = 3; % use 3 for single day, 3 hunting session experiments

if exp == 1
    ctl_mat = ctl_t2c_mat;
    dr_mat = dr_t2c_mat;

elseif exp == 2 
    ctl_mat = ctl_lat_mat;
    dr_mat = dr_lat_mat;
elseif exp == 3
    ctl_mat = ctl_ad_mat;
    dr_mat = dr_ad_mat;
end 

DLC_huntplot(ctl_mat(:,1:sess), dr_mat(:,1:sess), type, exp, hs)


