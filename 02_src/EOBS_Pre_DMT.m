
function [daily_pre, DMT, Date]=EOBS_Pre_DMT(points)

% this function is to get the daily precipitation for several grids or a
% region

% data from: https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php#datafiles
% Europe daily pre temperature, 0.25 gird
% 1950-01-01 -  2021-06-30
% eobs: 25N-71.5N x 25W-45E.
% Pre  

% updated on 2022.10.19
% new version of the EOBS data: v26 to 2022.6.31


path1='E:\Raw Data Series\EOBS\rr_ens_mean_0.25deg_reg_v26.0e.nc';
% ncdisp(path1)
lon=ncread( path1, 'longitude');
lat=ncread( path1, 'latitude');
time=ncread(path1, 'time');  

path2='E:\Raw Data Series\EOBS\tn_ens_mean_0.25deg_reg_v26.0e.nc'; % ncdisp(path2)
path3='E:\Raw Data Series\EOBS\tx_ens_mean_0.25deg_reg_v26.0e.nc';

N=size(points,1);
daily_pre=nan(N,length(time) ); daily_Tmax= nan(N,length(time) ); daily_Tmin= nan(N,length(time) );
if size(points, 2) ==2 % several grids
for n=1:N
    
    llon=points(n,1); llat=points(n,2);
    slon1=find(lon>llon-0.125 & lon<= llon+0.125) ;%经度
    slat1=find( lat>llat-0.125 & lat<=llat+0.125); %注意边界
    %     lon(slon1)
    %     lat(slat1)
    %load
    daily_pre(n,:)= squeeze( ncread( path1, 'rr', [slon1,slat1,1], [1,1, length(time)] ) )  ; % longitude,latitude,time, mm
    daily_Tmax(n,:)= squeeze( ncread( path3, 'tx', [slon1,slat1,1], [1,1, length(time)] ) ); % longitude,latitude,time, minimum temperature, degree
    daily_Tmin(n,:)= squeeze( ncread( path2, 'tn', [slon1,slat1,1], [1,1, length(time)] ) ); % longitude,latitude,time, minimum temperature, degree
end
% for output
t =datetime(1950,01,01)+caldays(0: size(time,1)- 1 )';
Date=[year(t), month(t), day(t) ];
daily_pre=[Date, daily_pre'];
DMT=[Date, mean( [daily_Tmax; daily_Tmin] )' ] ;

elseif  size(points, 2) ==4 % a region range: [lon1 left, lon2 right, lat south, lat north]
   
    % llon=points(n,1); llat=points(n,2);
    slon1= find(lon>points(1)-0.125 & lon<= points(1) +0.125);%经度
    slon2 =  find(lon>points(2)-0.125 & lon<= points(2) +0.125);%经度
    num_slon=slon2-slon1 +1;
    
    slat1=find( lat>points(3)-0.125 & lat<=points(3)+0.125); %注意边界
    slat2=find( lat>points(4)-0.125 & lat<=points(4)+0.125); %注意边界
    num_slat=slat2-slat1 +1;
    %     lon(slon1)
    %     lat(slat2) 
    %load
        daily_pre= squeeze( ncread( path1, 'rr', [slon1,slat1,1], [num_slon,num_slat, length(time)] ) ); % longitude,latitude,time, mm
    daily_Tmax= squeeze( ncread( path3, 'tx', [slon1,slat1,1], [num_slon,num_slat, length(time)] ) ); % longitude,latitude,time, minimum temperature, degree
    daily_Tmin= squeeze( ncread( path2, 'tn', [slon1,slat1,1], [num_slon,num_slat, length(time)] ) ); % longitude,latitude,time, minimum temperature, degree
    DMT=(daily_Tmax+daily_Tmin)/2;
    clear daily_Tmax daily_Tmin
    
end
    
% for output
t =datetime(1950,01,01)+caldays(0: size(time,1)- 1 )';
Date=[year(t), month(t), day(t) ];
