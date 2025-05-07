ephys_beh_dir = pwd; % run pipeline in Ephys Behavior directory containing code and the folder "animal_data"
addpath(genpath(ephys_beh_dir))

%% Ephys Rat Behavior Calculations

rat_times = dir('animal_data/*.csv');
ephys_anims = struct2cell(rat_times);

%Hunt Calculations
[t2c_mat,  lat_mat,  ad_mat] = DLC_huntcalc(ephys_anims, rat_times);

%% Ephys Rat Behavior Visualization
cd(ephys_beh_dir)
type = 2; % use 1 for day averages, use 2 for session averages
exp = 1; % 1 = time to capture, 2 = latency to attack, 3 = pursuit duration
hs = 4; % use 0 if for no plotting, or 1-3 for those sessions alone, or 4 for all
sess = 3; % use 3 for single day, 3 hunting session experiments

if exp == 1
    anim_mat = t2c_mat;
elseif exp == 2 
    anim_mat = lat_mat;
elseif exp == 3
    anim_mat = ad_mat;
end 

DLC_huntplot_1c(anim_mat(:,1:sess), type, exp, hs)

%% Analysis of Sex Differences
% Ephys Rat Behavior Calculations
% Males
m_times = dir('animal_data/Males/*.csv'); 
m_anims = struct2cell(m_times);

% Females
f_times = dir('animal_data/Females/*.csv'); 
f_anims = struct2cell(f_times);

%Hunt Calculations
[m_t2c_mat,  m_lat_mat,  m_ad_mat] = DLC_huntcalc(m_anims, m_times);
[f_t2c_mat,  f_lat_mat,  f_ad_mat] = DLC_huntcalc(f_anims, f_times);

if exp == 1
    m_mat = m_t2c_mat;
    f_mat = f_t2c_mat;

elseif exp == 2 
    m_mat = m_lat_mat;
    f_mat = f_lat_mat;
elseif exp == 3
    m_mat = m_ad_mat;
    f_mat = f_ad_mat;
end 

DLC_huntplot_mf(m_mat(:,1:sess), f_mat(:,1:sess), type, exp, hs)