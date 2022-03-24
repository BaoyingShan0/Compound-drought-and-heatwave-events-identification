% The aims of this demo include:
% (1) Calculating daily SPI and SHI value (non stationarity considered)
% (2) Getting the threholds for removing and merging with the goals that drought and heatwave events should extreme and independent. 
% (3) Identifing compound drought and heatwave events
% (4) Plotting figures

% What data you need is the daily precipitation(mm), daily max temperature(¡æ), daily main temperature(¡æ);

% This demo is written by Baoying Shan (Ghent University), 2022.3.24

%% Load data
clc; clear; close all;
load data.mat; % Date=[year, month, day], Pre=[daily pre], DMT=[daily mean temperature averaged of maximum and minimum]

% % where the Tmin_Tmax_Tave_ET0_pre data from? add the path
T=Tmin_Tmax_Tave_ET0_pre( :,1:5 );
Pre=Tmin_Tmax_Tave_ET0_pre( :,8 );
Date=Tmin_Tmax_Tave_ET0_pre( :, 1:3); % [year, month, day]

%%  366 days ---> 365 days
mm = find( Date(:,2)==2 & Date(:,3)==29);
Pre(mm-1)=Pre(mm-1) + Pre(mm); % for pre, sum is better
T(mm-1, 4:5)=( T(mm-1, 4:5)+T(mm, 4:5) )/2; % for temperature, sum is better
Date(mm,:)=[]; Pre(mm,:)=[];  T(mm,:)=[];

DMT=( T(:,4)+T(:,5))/2; % Daily mean temperature: the average of daily maximum and minimum temperatures
N=size(Date,1);
t1=datetime(Date);
years=Date(end,1)-Date(1,1)+1;


%% Calculation of SPI (Standardized Precipitation Index)
scale_p=15; % GIven the cummulation period for drought
% Assum the precipitation time series is stationary
[  spi_sta, p_h, p_sim_cdf , p_emp_cdf, p_NSE , p_best]= SPI_best_dist ( Date, Pre, scale_p, "" ) ; % pick the best dis to fit pre and get SPI based on AIC and ks test, include the details about how to handle 0
sum( spi_sta==-inf | spi_sta==inf ) % check there is no infinite value of SPI
sum( p_h ) % check picked best distribution all past ks test

%%   Calculation of  SHI (Standardized Heatwave Index)
tic
scale_T=3; % GIven the cummulation period for heatwave
[shi, h, sim_cdf, emp_cdf , NSE, best]=SHI_d_best(Date, DMT, scale_T, "", 30);  % 30 is the 30 years window because of the nonstationary of temperature, the normal condition is the past 30 years.
sum( shi==-inf | shi==inf )
toc % could spend longer time, about 1000 seconds

SPI=spi_sta ;
SHI = shi;
DATE = Date;
years = 120;
t = datetime(DATE);

%% Plot  SPI and SHI
year_range=[1901, 1905];
time_lim=[t( 365*( year_range(1) -DATE(1,1) )+1), t( 365*( year_range(2) -DATE(1,1)) +365) ] ;
figure(1); ax(1)=subplot(2,1,1);  plot(t, SPI );
hold on; plot(t, -1*ones(length(t), 1 ), 'r' );  plot(t, 0*ones(length(t), 1 ), 'k' ); xlim( time_lim ); datetick('x','yyyy-mm' );  hold off;
legend('Drought index', '-1'); ylabel(' SPI '); grid on;

ax(2)=subplot(2,1,2);
ba=plot(t, SHI );
hold on; plot(t, 1*ones(length(t), 1 ), 'k' ); hold off;
xlim( time_lim ); ylim([-3,3]); datetick('x','yyyy-mm');
legend('Heatwave index', '1' ); ylabel(' SHI '); grid on;
linkaxes(ax, 'x');

%% Get the thresholds for removing and merging

start_th_d=-1; end_th_d=-1; start_th_h=1;  end_th_h=1; % Given the start and end threholds

% get a rought boundary for initial thresholds: the annual days in drought
% should less then 120 days and total drought number should large than 10
boundarys_d = bound_remove_merge("Dr", DATE, SPI, start_th_d, end_th_d,  [0: 1: 4 * scale_p],  120, 10  );
REMO_d= [boundarys_d(1):1:boundarys_d(2)]; MERG_d = [-inf, 0:1:boundarys_d(2) ];

boundarys_h = bound_remove_merge("Hw", DATE, SHI, start_th_h, end_th_h,  [0: 1: 4 * scale_p],  120, 10  );
REMO_h= [boundarys_h(1):1:boundarys_h(2)]; MERG_h = [-inf, 0:1:boundarys_h(2) ];

p=0.05; % to make sure the best distribution pass the ks test at 0.05 significant level

[ best_com_dr , table_dr]= remo_merg_dr( DATE, SPI, REMO_d, MERG_d,   p, start_th_d, end_th_d ) ;
table_dr.("Scale") =  ones ( height(table_dr), 1) * scale_p;
table_dr.("Start_end_thresholds") =  ones ( height(table_dr), 1) * start_th_d;
best_com_dr


[ best_com_hw , table_hw]=remo_merg_hw( DATE, SHI , REMO_h, MERG_h, p, start_th_h, end_th_h) ;
table_hw.("Scale") =  ones ( height(table_hw), 1) * scale_T;
table_hw.("Start_end_thresholds") =  ones ( height(table_hw), 1) * start_th_h;
best_com_hw

%% Drought and heatwave events

thresholds =[ str2double( best_com_dr(1,1:2) );str2double( best_com_hw(1,1:2) ) ] ;

drought_daily=PRM2_drought_identification( DATE, SPI, start_th_d,end_th_d, thresholds(1,1), thresholds(1,2) );
% [year, month,  day, SPI, columns 5 to 6 are results after pre
% identification, columns 7 to 8 are results after removing, columns 9 to 10 are results after merging];
% 1 in columns 5, 7, 9 mean this day in drought; negative number in columns
% 6, 8, 10 are  the drought numbers.
heatwave_daily=PRM2_heatwave_identification( DATE, SHI, start_th_h, end_th_h, thresholds(2,1), thresholds( 2,2) );

size_dr=size(drought_daily, 2); size_hw=size(heatwave_daily, 2);
SPI_0=SPI;  SPI_0(drought_daily(:, size_dr-1)==0 )=nan;  % non-dry days are nan
SHI_0=SHI;  SHI_0(heatwave_daily(:, size_hw-1)==0 )=nan;

% daily to events
dr = daily_2_events(drought_daily); % [Number, duration, severity, magnitude,  the duration of the non-dry period, the inter-arrival time]
hw = daily_2_events(heatwave_daily, 'hw');

% annual nummber of events and annual days in events
number_of_dr=size(dr,1)/years;
annual_days_dr=sum( dr(:,2) )/years;
number_of_hw=size(hw,1)/years;
annual_hw=sum( hw(:,2) )/years;

%   hist for severity, duration, intensity
figure(1); % hist for severity, duration, intensity
subplot(2, 3,1); histogram(dr(:,2)); title('Histogram of duration for drought');  grid on;
subplot(2, 3,2); histogram( dr(:,3) ); title('Histogram of severity for drought'); grid on;
subplot(2, 3,3); histogram(dr(:,4)); title('Histogram of intensity for drought'); grid on;
subplot(2, 3,4); histogram(hw(:,2)); title('Histogram of duration for heatwave');  grid on;
subplot(2, 3,5); histogram(hw(:,3)); title('Histogram of severity for heatwave'); grid on;
subplot(2, 3,6); histogram(hw(:,4)); title('Histogram of intensity for heatwave'); grid on;

% scatter plot between duration-severity, duration-intensity, severity-intensity
figure(2);
subplot(2, 3,1); plot(dr(:,2), dr(:, 3) , '.');  xlabel('Duration (days)'); ylabel('Severity'); title('Dr');
subplot(2, 3,2); plot(dr(:,2), dr(:, 4) , '.');  xlabel('Duration (days)'); ylabel('Intensity'); title('Dr');
subplot(2, 3,3); plot(dr(:,4), dr(:, 3) , '.');  ylabel('Severity'); xlabel('Intensity'); title('Dr');
subplot(2, 3,4); plot(hw(:,2), hw(:, 3) , '.');  xlabel('Duration (days)'); ylabel('Severity'); title('Hw');
subplot(2, 3,5); plot(hw(:,2), hw(:, 4) , '.');  xlabel('Duration (days)'); ylabel('Intensity'); title('Hw');
subplot(2, 3,6); plot(hw(:,4), hw(:, 3) , '.');  ylabel('Severity'); xlabel('Intensity'); title('Hw');


%% Plot  drought and heatwave events
year_range=[1975, 1976];
time_lim=[t( 365*( year_range(1) -DATE(1,1) )+1), t( 365*( year_range(2) -DATE(1,1)) +365) ] ;

%time_lim=[datetime(2018,06,01), datetime(2018,07,31)]
figure(2);
ax(1)=subplot(2,1,1);  plot(t, drought_daily(:,4) );  hold on;
cc=drought_daily(:, 4); cc(drought_daily(:, end-1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', [1, 0.8, 0 ]  );  ba(1).BarWidth=1;
plot(t, start_th_d*ones(length(t), 1 ), 'r--' ); ylim([-3,3])
title( ['threshold for removing = '  num2str(thresholds(1,1))  ' days, threshold for merging =  ' num2str(thresholds(1,2)) ])
xlim( time_lim ); datetick('x','yyyy-mm' );   hold off;  ylabel(' SPI '); grid on; % legend('Drought index', 'Start threshhold');

ax(2)=subplot(2,1,2);
plot(t, heatwave_daily(:,4), 'k' , 'LineWidth', 0.8);   hold on;
cc=heatwave_daily(:, 4); cc(heatwave_daily(:, end-1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', [0.6, 0.8, 0.5 ]  );   ba(1).BarWidth=1;
plot(t, start_th_h*ones(length(t), 1 ), 'r--' ); hold off; ylabel(' SHI '); grid on;
title( ['threshold for removing = '  num2str(thresholds(2,1))  ' days, threshold for merging =  ' num2str(thresholds(2,2)) ])
xlim( time_lim );  ylim( [-3, 3] ); set(gca, 'YTick', -3:1:3); datetick('x','yyyy-mm');

linkaxes(ax, 'x');

%% Identification of CDHW
for i=1:4
    type=i;
    [compound, compound_events]=identify_comound_all(SPI_0, SHI_0, SPI, SHI, type); 
    % compound events: [duration, marginal severity of drought, marginal severity of heatwave]
    
end

%% Plot drought period, heatwave period, and compound period

types=[  "Union"    "Hw|Dr"    "Dr|Hw"    "Intersection"];
year_range=[2001, 2020];
time_lim=[t( 365*( year_range(1) -DATE(1,1) )+1), t( 365*( year_range(2) - DATE(1,1) ) +365) ] ;

figure;
ax(1)=subplot(6,1,1);  plot(t, SPI ,  'color', [0.85,0.33,0.10]);  hold on;
cc=drought_daily(:, 4); cc(drought_daily(:, end - 1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', [1, 0.8, 0 ]  );  ba(1).BarWidth=1;
plot(t, start_th_d*ones(length(t), 1 ), 'r' );
xlim( time_lim ); datetick('x','yyyy-mm' );   hold off;  ylabel(' SPI '); grid on; title('Drought');
% legend('SPI', 'Drought duration', '-1');

ax(2)=subplot(6,1,2);
hold on;
cc=heatwave_daily(:, 4); cc(heatwave_daily(:, end - 1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', 'r');   ba(1).BarWidth=1; %, [0.6, 0.8, 0.5 ]
plot(t, SHI,  'color' ,[0.47,0.5,0.5] );
plot(t, start_th_h*ones(length(t), 1 ), 'r' ); hold off;
xlim( time_lim );  ylim( [-3, 3] ); datetick('x','yyyy-mm'); title('Heatwave')
ylabel(' SHI '); grid on;
% legend('SHI', 'Heatwave duration', '1');

for i=1:4
    type=i;
    [compound, compound_events]=identify_comound_all(SPI_0, SHI_0, SPI, SHI, type); 
    ax(i+2)=subplot(6,1,i+2); % for compound period
    ba= bar(t, compound(:,1));   ba(1).BarWidth=2;
    xlim( time_lim ); datetick('x','yyyy-mm');  yticks([0, 1]);
    title(types(type) );
    %legend('Identified compound duration')
end
xlabel('Date')
linkaxes(ax, 'x');

% Results are save to results.mat