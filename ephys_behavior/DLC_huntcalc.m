function [avg_t2c_mat, avg_lat_mat, avg_ad_mat] =  DLC_huntcalc(cohort, directory)


%% Make Structs
avg_t2c = {};
avg_lat = {};
avg_attack_dur = {};
anim_sum = {};

fps = 20;

%% Main Calculation Loop
for i = 1:size(cohort,2)
 
t2c = [];
lat = [];
attack_dur = [];
    
%Grab Animal Data
rat_mat = readmatrix(cohort{1,i});
anim_num = char(extractBetween(directory(i).name,'_','.'));

% Variables
n_sessions = 1:max(rat_mat(:,3));
n_days = max(rat_mat(:,2));


%Total Time to Capture
t2c(:,1) = rat_mat(:,3);
t2c(:,2) = (rat_mat(:,7) - rat_mat(:,5))./fps;

%Latency
lat(:,1) = rat_mat(:,3);
lat(:,2) = (rat_mat(:,6) - rat_mat(:,5))./fps;

%Attack Duration
attack_dur(:,1) = rat_mat(:,3);
attack_dur(:,2) = (rat_mat(:,7) - rat_mat(:,6))./fps;

%Parse into HS
max_sessions = 9; 

HS_lat = NaN(13,max_sessions);
HS_attackdur = NaN(13,max_sessions);
HS_t2c = NaN(13,max_sessions);

    for ii = 1:max_sessions,
        
        crik = find(rat_mat(:,3) == ii);
        if isempty(crik)
            continue
        else
        crik1 = crik(1);
        criklast = crik(end);
        n_val = length(lat(crik1:criklast,2));

        HS_lat(1:n_val,ii) = lat(crik1:criklast,2);
        HS_attackdur(1:n_val,ii) = attack_dur(crik1:criklast,2);
        HS_t2c(1:n_val,ii) = t2c(crik1:criklast,2);
        end
    end

 % Averages
    avg_lat{1,i} = anim_num;
    avg_lat{2,i} = nanmean(HS_lat,1);
    
    avg_attack_dur{1,i} = anim_num;
    avg_attack_dur{2,i} = nanmean(HS_attackdur,1);
    
    avg_t2c{1,i} = anim_num;
    avg_t2c{2,i} = nanmean(HS_t2c,1);
    
  % Summary Cell
    anim_sum{1,i} = anim_num;
    anim_sum{2,i} = HS_t2c;
    anim_sum{3,i} = HS_lat;
    anim_sum{4,i} = HS_attackdur;
end

%% Extract Data for Averages
n_animals = size(cohort,2);
avg_lat_mat = NaN(n_animals,max_sessions);
avg_ad_mat = NaN(n_animals,max_sessions);
avg_t2c_mat = NaN(n_animals,max_sessions);

for iii = 1:n_animals
    
    if length(avg_lat{2,iii}) == 3
        avg_lat_mat(iii,1:3) = avg_lat{2,iii};
        avg_ad_mat(iii,1:3) = avg_attack_dur{2,iii};
        avg_t2c_mat(iii,1:3) = avg_t2c{2,iii};
    else
        avg_lat_mat(iii,:) = avg_lat{2,iii};
        avg_ad_mat(iii,:) = avg_attack_dur{2,iii};
        avg_t2c_mat(iii,:) = avg_t2c{2,iii};
        
    end
end

%% Average by HS


hsavg_lat = mean(avg_lat_mat,1);
hsavg_ad = mean(avg_ad_mat,1);
hsavg_t2c = mean(avg_t2c_mat,1);


%% Average by Day
d1avg_t2c = mean(avg_t2c_mat(:,1:3),2);

d2avg_t2c = mean(avg_t2c_mat(:,4:5),2);
d3avg_t2c = mean(avg_t2c_mat(:,6:7),2);
d4avg_t2c = mean(avg_t2c_mat(:,8:9),2);
davg_t2c = [d1avg_t2c d2avg_t2c];

d1avg_ad = mean(avg_ad_mat(:,1:3),2);
d2avg_ad = mean(avg_ad_mat(:,4:5),2);
d3avg_ad = mean(avg_ad_mat(:,6:7),2);
d4avg_ad = mean(avg_ad_mat(:,8:9),2);
davg_ad = [d1avg_ad d2avg_ad];

d1avg_lat = mean(avg_lat_mat(:,1:3),2);
d2avg_lat = mean(avg_lat_mat(:,4:5),2);
d3avg_lat = mean(avg_lat_mat(:,6:7),2);
d4avg_lat = mean(avg_lat_mat(:,8:9),2);
davg_lat = [d1avg_lat d2avg_lat];

end