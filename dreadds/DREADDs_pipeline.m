%% DREADDs PIPELINE for Leman et al., 2025

dr_code_dir = pwd;% uigetdir(); % Choose directory containing DREADDs code
addpath(genpath(DR_code_dir))

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
type = 2;
exp = 1;
hs = 4;
sess = 3;

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


