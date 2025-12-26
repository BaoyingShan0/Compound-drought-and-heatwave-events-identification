
%% Introduction
% This script implements a proposed statistical method for identifying drought and heatwave events from climate data. 
% It includes data pre-processing, calculation of SPI and SHI, events identification, and visualization.

% Please cite this paper if you use this method:  
% Shan, B., Verhoest, N. E. C., and De Baets, B.: Identification of compound drought and heatwave events on a daily scale and across four seasons, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2023-147, 2023.

% Auther: Baoying Shan, Gent University, (baoying.shan@ugent.be)
% Date: 2024.03.24


% This script includes
% part 1: Load functions and data
% part 2: Calculate SPI and SHI
% Part 3: Removal and merging processes
% Part 4: Identify drought and heatwave events
% Part 5: Identify compound drought and heatwave events
% Part 6: Plotting

%% Part 1: Load functions and data

addpath(genpath('./02_src'));
load('01_data/data0.mat'); 
% which saved Pre(the daily precipitation), DMT(daily mean temperature), and Date (corresponding year, month, day) 

% Data pre-processing:
% Convert 366 days to 365 days
mm = find(Date(:,2) == 2 & Date(:,3) == 29);
Pre = pre; Pre(mm - 1, :) = pre(mm - 1, :) + pre(mm, :);
DMT = dmt; DMT(mm - 1, :) = (dmt(mm - 1, :) + dmt(mm, :)) / 2;
Date(mm, :) = []; Pre(mm, :) = []; DMT(mm, :) = [];
t = datetime(Date);
years = Date(end, 1) - Date(1, 1) + 1;

%% Part 2: compariring Calculate SPI and SHI

NSP = 30; % parameter for non-stationarity
scale_p = 30; % accumulation period for droughts
scale_T = 3; % accumulation period for heatwaves

% SPI nonparametric
SPI_nonpar = SI_nonparametric(Date, Pre, scale_p); % 0.26 seconds
% SPI best distribution 
tic; SPI_best_dist = SPI_best_dist(Date, Pre, scale_p, NSP); toc; % 34 seconds
plot(SPI_nonpar, SPI_best_dist,'.'); xlabel('SPI_nonpar');ylabel("SPI_best_dist");
grid on; hold on; plot(-3:0.1:3, -3:0.1:3,'r--');


%  SHI
tic
[SHI_nonpar] = SI_nonparametric(Date, DMT, scale_T);
toc
tic 
%  SHI
[SHI_best_dist] = SHI_best_dist(Date, DMT, scale_T, NSP);
toc
figure; scatter(SHI_nonpar, SHI_best_dist); xlabel('SHI_nonpar');ylabel("SHI_best_dist");
grid on; hold on; plot(-3:0.1:3, -3:0.1:3,'r--');
%% Part 3: Removal and merging processes
% Parameter setting
start_th_d=-1; end_th_d=-1; % pre-identification thresholds for droughts
start_th_h=1;  end_th_h=1; %  pre-identification thresholds for heatwaves
p=0.05;
events_number_control=0.1; events_days_control=120;


% Get removal and merging thresholds for droughts
tic
drought_removal_merging_best_dist =  remo_merg("d",Date, SPI_best_dist, ...
    scale_p, p, start_th_d, end_th_d, events_number_control, events_days_control); % "d" represents droughts
drought_removal_threshold_best_dist = cc(1)  
drought_merging_threshold_best_dist = cc(2)
toc

tic
drought_removal_merging_nonpar =  remo_merg("d",Date, SPI_nonpar, ...
    scale_p, p, start_th_d, end_th_d, events_number_control, events_days_control); % "d" represents droughts
drought_removal_threshold_nonpar= drought_removal_merging_nonpar(1)  
drought_merging_threshold_nonpar = drought_removal_merging_nonpar(2)
toc

% Get removal and merging thresholds for heawtaves
heatwave_removal_merging_best_dist = remo_merg("h",Date, SHI_best_dist, ...
    scale_T, p, start_th_h, end_th_h, events_number_control, events_days_control); % "h" represnets heatwaves
heatwave_removal_threshold_best_dist=heatwave_removal_merging_best_dist(1) 
heatwave_merging_threshold_best_dist=heatwave_removal_merging_best_dist(2)

heatwave_removal_merging_nonpar = remo_merg("h",Date, SHI_nonpar, ...
    scale_T, p, start_th_h, end_th_h, events_number_control, events_days_control); % "h" represnets heatwaves
heatwave_removal_threshold_nonpar=heatwave_removal_merging_nonpar(1); 
heatwave_merging_threshold_nonpar=heatwave_removal_merging_nonpar(2);
heatwave_removal_merging_nonpar

%% Part 4: Identify drought and heatwave events
drought_daily_best_dist = PRM_extreme_identification("d", Date, SPI_best_dist,start_th_d,end_th_d,...
    drought_removal_threshold_best_dist, drought_merging_threshold_best_dist ); % the ninth column records if this day is in a drought (1 for yes, 0 for no)
drought_events_best_dist = daily_2_events(drought_daily_best_dist, "d"); % drought events


drought_daily_nonpar = PRM_extreme_identification("d", Date, SPI_nonpar,start_th_d,end_th_d,...
    drought_removal_threshold_nonpar, drought_merging_threshold_nonpar ); % the ninth column records if this day is in a drought (1 for yes, 0 for no)
drought_events_nonpar = daily_2_events(drought_daily_nonpar, "d"); % drought events

% taking year 1976 as an example
year_range=[1950, 1955];
time_lim=[t( 365*( year_range(1) -Date(1,1) )+1), t( 365*( year_range(2) - Date(1,1) ) + 365) ] ;
figure(1);
% plot SPI and drought identification results
ax(1)=subplot(2,1, 1);
hold on;
cc=drought_daily_best_dist(:, 4); cc(drought_daily_best_dist(:, end - 1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', [1, 0.8, 0 ]  );  ba(1).BarWidth=1;
plot(t, SPI_best_dist, '.-',  'color', [0.85,0.33,0.10], 'LineWidth', 0.8 );
plot(t, start_th_d*ones(length(t), 1 ), 'r--',  'LineWidth', 0.5);
xlim( time_lim );hold off; ylabel(' SPI '); grid on; ti=title('Drought'); set(ti , 'FontWeight', 'normal' )
datetick('x','yyyy-mm' );   set(gca, 'xticklabel', []); box on;

ax(2)=subplot(2,1, 2);
hold on;
cc=drought_daily_nonpar(:, 4); cc(drought_daily_nonpar(:, end - 1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', [1, 0.8, 0 ]  );  ba(1).BarWidth=1;
plot(t, SPI_nonpar, '.-',  'color', [0.85,0.33,0.10], 'LineWidth', 0.8 );
plot(t, start_th_d*ones(length(t), 1 ), 'r--',  'LineWidth', 0.5);
xlim( time_lim );hold off; ylabel(' SPI '); grid on; ti=title('Drought'); set(ti , 'FontWeight', 'normal' )
datetick('x','yyyy-mm' );   set(gca, 'xticklabel', []); box on;
linkaxes(ax, 'x');

%% heatwave
heatwave_daily_best_dist = PRM_extreme_identification("h", Date, SHI_best_dist,start_th_h,end_th_h,...
    heatwave_removal_threshold_best_dist, heatwave_merging_threshold_best_dist ); % the ninth column records if this day is in a heatwave (1 for yes, 0 for no)
heatwave_events_best_dist = daily_2_events(heatwave_daily_best_dist, "h"); % heatwave events


heatwave_daily_nonpar = PRM_extreme_identification("h", Date, SHI_nonpar,start_th_h,end_th_h,...
    heatwave_removal_threshold_nonpar, heatwave_merging_threshold_nonpar ); % the ninth column records if this day is in a heatwave (1 for yes, 0 for no)
heatwave_events_nonpar = daily_2_events(heatwave_daily_nonpar, "h"); % heatwave events

% taking year 1976 as an example
year_range=[1950, 1955];
time_lim=[t( 365*( year_range(1) -Date(1,1) )+1), t( 365*( year_range(2) - Date(1,1) ) + 365) ] ;
figure(1);
% plot SHI and heatwave identification results
ax(1)=subplot(2,1, 1);
hold on;
cc=heatwave_daily_best_dist(:, 4); cc(heatwave_daily_best_dist(:, end - 1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', [1, 0.8, 0 ]  );  ba(1).BarWidth=1;
plot(t, SHI_best_dist, '.-',  'color', [0.85,0.33,0.10], 'LineWidth', 0.8 );
plot(t, start_th_h*ones(length(t), 1 ), 'r--',  'LineWidth', 0.5);
xlim( time_lim );hold off; ylabel(' SHI '); grid on; ti=title('heatwave'); set(ti , 'FontWeight', 'normal' )
datetick('x','yyyy-mm' );   set(gca, 'xticklabel', []); box on;

ax(2)=subplot(2,1, 2);
hold on;
cc=heatwave_daily_nonpar(:, 4); cc(heatwave_daily_nonpar(:, end - 1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', [1, 0.8, 0 ]  );  ba(1).BarWidth=1;
plot(t, SHI_nonpar, '.-',  'color', [0.85,0.33,0.10], 'LineWidth', 0.8 );
plot(t, start_th_h*ones(length(t), 1 ), 'r--',  'LineWidth', 0.5);
xlim( time_lim );hold off; ylabel(' SHI '); grid on; ti=title('heatwave'); set(ti , 'FontWeight', 'normal' )
datetick('x','yyyy-mm' );   set(gca, 'xticklabel', []); box on;
linkaxes(ax, 'x');


%% Part 5: Identify compound events

types=["d-and-h", "d-or-h", "d-cond-h", "h-cond-d"];
[~, compound_daily]=identify_compound( drought_daily, heatwave_daily, 1); % 1 represents "d-and-h"

%% Part 6: Plotting

% taking year 1976 as an example
year_range=[1976, 1976];
time_lim=[t( 365*( year_range(1) -Date(1,1) )+1), t( 365*( year_range(2) - Date(1,1) ) + 365) ] ;

figure(1);

% plot SPI and drought identification results
ax(1)=subplot(8,1, [1:2]);
hold on;
cc=drought_daily(:, 4); cc(drought_daily(:, end - 1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', [1, 0.8, 0 ]  );  ba(1).BarWidth=1;
plot(t, SPI_best_dist, '.-',  'color', [0.85,0.33,0.10], 'LineWidth', 0.8 );
plot(t, start_th_d*ones(length(t), 1 ), 'r--',  'LineWidth', 0.5);
xlim( time_lim );hold off; ylabel(' SPI '); grid on; ti=title('Drought'); set(ti , 'FontWeight', 'normal' )
datetick('x','yyyy-mm' );   set(gca, 'xticklabel', []); box on;

% plot SHI and heatwave identification results
ax(2)=subplot(8,1, [3:4] );
hold on;
cc=heatwave_daily(:, 4); cc(heatwave_daily(:, end - 1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', [0.6, 0.8, 0.5 ]  );   ba(1).BarWidth=1;
plot(t, SHI_best_dist, '.-', 'color', [0.47,0.67,0.19],   'LineWidth', 0.8);
plot(t, start_th_h*ones(length(t), 1 ), 'r--',  'LineWidth', 0.5 ); hold off;
xlim( time_lim );  ylim( [-3, 3] );  ti=title('Heatwave'); set(ti , 'FontWeight', 'normal' )
ylabel(' SHI '); grid on; datetick('x','yyyy-mm');  set(gca, 'xticklabel', []); box on;

% Compound events identification results
for i=1:4
    type=i;
    compound=identify_compound(drought_daily, heatwave_daily, type);

    ax=subplot( 8, 1, 4+i ); % for compound period
    aa=get(ax,'position');
    ax=subplot(8,1,4+i , 'Position',aa+[0,0,0,0]);
    ba= bar(t, compound(:,1));   ba(1).BarWidth=2; ba(1).FaceColor=[0.5,0.5,0.5];
    xlim( time_lim );  yticks([0, 1]); ylabel('CDHW'); set(gca, 'yticklabel', [])
    ti = title(types(type)); set(ti , 'FontWeight', 'normal' );grid on;

    if type<4
        datetick('x','yyyy-mm');
        set(gca, 'xticklabel', []);
    else
        datetick('x','yyyy-mm'); xlabel('Date')
    end
end
linkaxes(ax, 'x');
set(gcf, 'Position', [150, 120, 1000, 700]);
% output
% print(gcf,'identification_results.png','-dpng','-r600');
