% 2022.10.19
% Create by Baoying Shan, baoying.shan@ugent.be

% An example to identify drought, heatwave and compound drought and
% heatwave events

% taking Be_Bra as an example, 
% 51��18'27.4"N, 4��31'11.4"E, BRASSCHAAT,  https://www.icos-belgium.be/ESBrasschaat.php 

clc;clear;
lon=4+31/60+11.4/3600; lat=51+18/60+27.4/3600;
points= [ lon, lat]; % latitude and longitude of the given site

%% load meteorological data
tic
[Pre0, DMT0]=EOBS_Pre_DMT(points);
toc

t =[datetime(1950,01,01): 1: datetime(2022,12,31)]';
Date=[year(t), month(t), day(t) ];
Pre= [ Pre0(:,4); nan( size(Date,1)- size(Pre0,1) , 1) ] ; % 1970-2022
DMT= [ DMT0(:,4);  nan(size(Date,1)- size(Pre0,1),1 ) ];

%%  366 ---> 365 days
mm = find( Date(:,2)==2 & Date(:,3)==29);
Pre( mm-1)=Pre( mm-1) + Pre( mm); % sum or average? for pre, maybe sum is better
DMT( mm-1, :)=( DMT( mm-1, :)+DMT( mm, :) )/2;
Date( mm,:)=[]; Pre(mm)=[];  DMT(mm)=[];

years = Date(end,1)-Date(1,1)+1;
N=size(Date,1);
N1=365*(30-1)+1;
t = datetime(Date);

%% SPI
scale_p=15%30; % cummulation period

% stationarity, 30s
tic
[spi_sta, p_h, p_sim_cdf , p_emp_cdf, p_NSE , p_best]= SPI_best_dist ( Date, Pre, scale_p,"" ) ; % pick the best dis to fit pre and get SPI based on AIC and ks test, include the details about how to handle 0
toc

% % non stationarity, for comparison
% [ spi_nonsta, p_h_nonsta, p_sim_cdf_nonsta, p_emp_cdf_nonsta, p_NSE_nonsta, p_best_nonsta ]=SPI_best_dist( Date, Pre, scale_p, "", 30 );  % 30 is the 30 years because of the nonstationary of temperature, the normal condition is the past 30 years.
% quality control
sum( spi_sta==-inf | spi_sta==inf ) % make sure there is no 0
sum( p_h )
max( nanmin( reshape( spi_sta, 365,[] ) , [], 2 ) ) % whether it is less than -1 to detect the zero numbers
SPI=spi_sta ;

% plot( t1, isnan(spi_sta) )

%%  SHI , 2000s,
tic
scale_T=3;
[shi, h, sim_cdf, emp_cdf , NSE, best]=SHI_d_best(Date, DMT, scale_T, "", 30);  % 30 is the 30 years because of the nonstationary of temperature, the normal condition is the past 30 years.
sum( shi==-inf | shi==inf )
toc
SHI = shi;
% The first 29 years are included in the next analysis


%% plot  SPI and SHI
year_range=[2016, 2018];
time_lim=[t( 365*( year_range(1) -Date(1,1) )+1), t( 365*( year_range(2) -Date(1,1)) +365) ] ;
figure(1);
ax(1)=subplot(2,1,1);  plot(t, SPI );
hold on; plot(t, -1*ones(length(t), 1 ), 'r' );  plot(t, 0*ones(length(t), 1 ), 'k' ); xlim( time_lim ); datetick('x','yyyy-mm' );  hold off;
legend('Drought index', '-1'); ylabel(' SPI '); grid on;
% title('xx')

ax(2)=subplot(2,1,2);
ba=plot(t, SHI ); %ba.FaceColor=[.93 .69 .13];
hold on; plot(t, 1*ones(length(t), 1 ), 'r' );  plot(t, 0*ones(length(t), 1 ), 'k' );  hold off;
xlim( time_lim ); ylim([-3,3]); datetick('x','yyyy-mm');
legend('Heatwave index', '1' ); ylabel(' SHI '); grid on;
linkaxes(ax, 'x'); 

%% removing and merging
start_th_d=-1; end_th_d=-1; start_th_h=1;  end_th_h=1; % set pre_identification threshsold

p=0.05; % to make sure the best distribution pass the ks test at 0.05 significant level


% scale check
if max( nanmin( reshape(SPI, 365, []), [],2 )  ) > start_th_d || min( nanmax( reshape(SHI, 365, []), [],2 )  ) < start_th_h
    error("The scales and pre identification thresholds are not consistent"); % there would be days without any days in dr or hw
end
% REMO_d=[45:-1: 0]; MERG_d=[45:-1: 0]; REMO_h=[20:-1: 0]; MERG_h=[20:-1:0];
% get a rought boundary for initial thresholds the annual days in drought
% should less then 120 days and total drought number over 91 years should
% large than 10
boundarys_d = bound_remove_merge("Dr", Date, SPI, start_th_d, end_th_d,  [0: 1: 4 * scale_p],  120, 10  );
REMO_d= [boundarys_d(1):1:boundarys_d(2)]; MERG_d = [-inf, 0:1:boundarys_d(2) ];

boundarys_h = bound_remove_merge("Hw", Date, SHI, start_th_h, end_th_h,  [0: 1: 4 * scale_T],  120, 10  );
REMO_h= [boundarys_h(1):1:boundarys_h(2)]; 
REMO_h=[1:1:10]
MERG_h = [-inf, 0:1:boundarys_h(2) ];


%Removing and merging for drought
[ best_com_dr , table_dr]= remo_merg_dr( Date, SPI, REMO_d, MERG_d,   p, start_th_d, end_th_d ) ;
table_dr.("Scale") =  ones ( height(table_dr), 1) * scale_p;
table_dr.("Start_end_thresholds") =  ones ( height(table_dr), 1) * start_th_d;
best_com_dr


% Removing and merging
[ best_com_hw , table_hw]=remo_merg_hw( Date, SHI , REMO_h, MERG_h, p, start_th_h, end_th_h) ;
table_hw.("Scale") =  ones ( height(table_hw), 1) * scale_T;
table_hw.("Start_end_thresholds") =  ones ( height(table_hw), 1) * start_th_h;
best_com_hw

% Same removing and merging procedures could be used to identify extreme wet and coldwave
%  extreme precipitation
% start_th_ep=1;  end_th_ep=1;
% [ best_com_ep , table_ep]= remo_merg_hw( Date, SPI, REMO_d, MERG_d,   0.05, start_th_ep, end_th_ep ) ;
% table_ep.("Scale") =  ones ( height(table_ep), 1) * scale_p;
% table_ep.("Start_end_thresholds") =  ones ( height(table_ep), 1) * start_th_ep;
% best_com_ep
% 
% 
% % remove merge for coldwave
% start_th_cw=-1;  end_th_cw=-1; 
% [ best_com_cw , table_cw]=remo_merg_dr( Date, SHI , REMO_h, MERG_h, 0.05, start_th_cw, end_th_cw) ;
% table_cw.("Scale") =  ones ( height(table_cw), 1) * scale_T;
% table_cw.("Start_end_thresholds") =  ones ( height(table_cw), 1) * start_th_cw;
% best_com_cw


%% identification of extremes��
ths_dr =str2double( best_com_dr(1,1:2) ) ;
drought_daily=PRM2_drought_identification( Date, SPI, start_th_d, end_th_d, ths_dr(1), ths_dr(2) );

ths_hw =str2double( best_com_hw(1,1:2) ) ;
heatwave_daily=PRM2_heatwave_identification( Date, SHI, start_th_h, end_th_h, ths_hw(1), ths_hw(2) );

% ths_ep =str2double( best_com_ep(1,1:2) ) ;
% extremepre_daily=PRM2_heatwave_identification( DATE, SPI, start_th_ep,end_th_ep, ths_ep(1,1), ths_ep(1,2)  );%17, 41);
% 
% ths_cw =str2double( best_com_cw(1,1:2) ) ;
% coldwave_daily=PRM2_drought_identification( DATE, SHI, start_th_cw,start_th_cw, ths_cw(1,1), ths_cw(1,2)  );%17, 41);

%% plot  extremes
year_range=[2017 2020];
time_lim=[t( 365*( year_range(1) -Date(1,1) )+1), t( 365*( year_range(2) -Date(1,1)) +365) ] ;

%time_lim=[datetime(2018,06,01), datetime(2018,07,31)]
figure(2);
ax(1)=subplot(2,1,1);  
plot(t, drought_daily(:,4) );  hold on;
cc=drought_daily(:, 4); cc(drought_daily(:, end-1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', [1, 0.8, 0 ]  );  ba(1).BarWidth=1;
% cc=extremepre_daily(:, 4); cc(extremepre_daily(:, end-1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', [0.5, 0.8, 0 ]  );  ba(1).BarWidth=1;
plot(t, start_th_d*ones(length(t), 1 ), '--', 'color',  [1, 0.8, 0 ] ); %ylim([-3,3]);
% plot(t, start_th_ep*ones(length(t), 1 ), '--', 'color',[0.5, 0.8, 0 ] ); ylim([-3,3]);
% title(['Thresholds for drought: [', num2str(ths_dr(1)), ',', num2str(ths_dr(2)), ']; ', 'Thresholds for extreme precipitation: [', num2str(ths_ep(1)), ',', num2str(ths_ep(2)), ']' ])
xlim( time_lim ); datetick('x','yyyy-mm' );   hold off;  ylabel(' SPI '); grid on; % legend('Drought index', 'Start threshhold');

ax(2)=subplot(2,1,2);
plot(t, heatwave_daily(:,4), 'k' , 'LineWidth', 0.8);   hold on;
cc=heatwave_daily(:, 4); cc(heatwave_daily(:, end-1)==0 )=nan;   ba= bar(t, cc , 'FaceColor',  [0.5, 0.2, 0.2 ] );   ba(1).BarWidth=1;
% cc=coldwave_daily(:, 4); cc(coldwave_daily(:, end-1)==0 )=nan;   ba= bar(t, cc , 'FaceColor', [0.6, 0.8, 0.5 ]   );  ba(1).BarWidth=1;
plot(t, start_th_h*ones(length(t), 1 ), '--', 'color',  [0.5, 0.2, 0.2 ]); 
% plot(t, start_th_cw*ones(length(t), 1 ), '--', 'color', [0.6, 0.8, 0.5 ]);
hold off; ylabel(' SHI '); grid on;
% title(['Thresholds for heatwave: [', num2str(ths_hw(1)), ',', num2str(ths_hw(2)), ']; ', 'Thresholds for coldwave: [', num2str(ths_cw(1)), ',', num2str(ths_cw(2)), ']' ])
xlim( time_lim );  ylim( [-3, 3] ); set(gca, 'YTick', -3:1:3); datetick('x','yyyy-mm'); xlabel('Date');
% legend('SHI','Heatwave', 'Coldwave', '1', '-1')
linkaxes(ax, 'x');

%% Daily to events

% If there is no revoming and merging
% drought_daily=PRM2_drought_identification( DATE, SPI, start_th_d,end_th_d, -inf, -inf );
% heatwave_daily=PRM2_heatwave_identification( DATE, SHI, start_th_h, end_th_h, -inf, -inf);

dr=nan(max( drought_daily(:,end) ), 1) ;
for i=1: max( drought_daily(:,end) ) %
    aa=drought_daily( drought_daily(:, end )==-i,:) ;
    dr(i,1:4)=[i, size(aa,1), -sum(aa(:,4) ), -sum(aa(:,4) )/size(aa,1)] ; %num, duration, severity, duration
    if i<2
        dr(i,5:6)=nan;
    else
        bb=find(drought_daily(:, end )==-i  ) ;
        cc=find( drought_daily(:, end )==-i+1 ); % last events
        dr(i,5)=bb(1)-cc(end)-1; % the duration of the non-dry period i
        dr(i,6)=dr(i, 2)+dr(i, 5); % the inter-arrival time
    end, end

% hw
hw=nan(max( heatwave_daily(:, end) ) ,1) ;
for i=1: max( heatwave_daily(:, end) ) %
    aa=heatwave_daily( heatwave_daily(:, end )==-i,:) ;
    hw(i,1:4)=[i, size(aa,1), sum(aa(:,4) ), sum(aa(:,4) )/size(aa,1)] ; %num, duration, severity, duration
    if i<2
        hw(i,5:6)=nan;
    else
        bb=find(heatwave_daily(:, end )==-i  ) ;
        cc=find( heatwave_daily(:, end )==-i+1 ); % last events
        hw(i,5)=bb(1)-cc(end)-1; % the duration of the non-heatwave period i
        hw(i,6)=hw(i, 2)+hw(i, 5); % the inter-arrival time
    end, end

% annual nummber of events and annual days in events
number_of_dr=size(dr,1)/years
annual_days_dr=sum( dr(:,2) )/years
number_of_hw=size(hw,1)/years
annual_hw=sum( hw(:,2) )/years

%% the identification of CDHW events 
types=[  " \itD\cup\itH"    "\itH|\itD"    "\itD|\itH"    "\itD\cap\itH"]; % means "Union", "Hw|Dr", "Dr|Hw", and "intersection", respectively
type=1;% union

SPI_0=SPI;  SPI_0(drought_daily(:, end-1)==0 )=nan;  % non-dry days are nan
SHI_0=SHI;  SHI_0(heatwave_daily(:, end-1)==0 )=nan;
compound=identify_comound(SPI_0, SHI_0, SPI, SHI, type); 

%% plot drought period, heatwave period, and compound period

year_range=[2015, 2021];
time_lim=[t( 365*( year_range(1) -Date(1,1) )+1), t( 365*( year_range(2) - Date(1,1) ) +365) ] ;

% figure;
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
%     types=["Union", "Hw|Dr", "Dr|Hw", "intersection"];
    compound=identify_comound(SPI_0, SHI_0, SPI, SHI, type); %
    number_com(i) = max( compound(:,2 ) )
    ax(i+2)=subplot(6,1,i+2); % for compound period
    ba= bar(t, compound(:,1));   ba(1).BarWidth=2;
    xlim( time_lim ); datetick('x','yyyy-mm');  yticks([0, 1]);
    title(types(type) ); 
    
end
xlabel( ['Date (lon= ', num2str(lon), ', lat= ', num2str(lat) , ')'])
linkaxes(ax, 'x');
legend('Compound events')
