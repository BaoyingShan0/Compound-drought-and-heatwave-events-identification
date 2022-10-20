%% 2022.9.1
%  put the compound daily by 4 compound ways in a big matrix
% save compound_daily.mat 
% % put all grids in a big matrix, whether this in compound events or not
% by four compound ways
% Com_hw_dr(i,j, type, :) = com_hw_dr(:,1) ;  Com_hw_ep(i,j,type, :) = com_hw_ep(:,1);
% Com_cw_dr(i,j,type,:) = com_cw_dr(:,1);  Com_cw_ep(i,j,type,:) = com_cw_ep(:,1);

clc;clear;
load uk_spi_shi.mat
load('uk_removing_merging 2.mat'); % including the heatwave_daily, drought_daily
DATE = Date;
years = Date(end,1)-Date(1,1)+1;
t = datetime(DATE);
start_th_d = -1; end_th_d = -1; start_th_h = 1;  end_th_h = 1;
start_th_cw = -1;  end_th_cw = -1; start_th_ep = 1;  end_th_ep = 1;
LLON=length(lon); LLAT=length(lat); TT=size(Date,1);

%% put the compound daily in a big matrix
types=[1,2,3,4];
[Com_hw_dr,  Com_hw_ep, Com_cw_dr, Com_cw_ep] = deal(nan(LLON, LLAT, 4,  N));
tic
non_nan=0;
% for all grids
for i=1:LLON
    i
    for j=1:LLAT
        
        % corresponding SPI, SHI, thresholds of this station
        SPI = squeeze( spis(i, j, : ) ); SHI = squeeze( shis( i, j, : ) );
        if ~(all(isnan(SPI)) & all(isnan(SHI)) )
            for type=1:4
                non_nan=non_nan+1;
                DMT = squeeze( DMT_all(i, j, : ) ); Pre= squeeze( Pre_all(i, j, : ) );
                
                heatwave_daily=squeeze( Heatwave_daily (i,j,:, 1) ); heatwave_daily( heatwave_daily==0 )=nan;
                drought_daily= squeeze( Drought_daily (i,j,:, 1 ) ); drought_daily( drought_daily==0 )=nan;
                coldwave_daily=squeeze( Coldwave_daily (i,j,:, 1 ) ); coldwave_daily( coldwave_daily==0 )=nan;
                extremepre_daily= squeeze(Extremepre_daily (i,j,:,1 ) ); extremepre_daily( extremepre_daily==0 )=nan;
                
                % compound hw and dr for each grid
                com_hw_dr= identify_comound ( heatwave_daily, drought_daily,SHI, SPI, type );
                % compound hw and ep for each grid
                com_hw_ep= identify_comound (heatwave_daily, extremepre_daily, SHI, SPI, type );
                % compound cw and dr for each grid
                com_cw_dr = identify_comound ( coldwave_daily, drought_daily, SHI, SPI, type );
                % compound cw and ep for each grid
                com_cw_ep= identify_comound (coldwave_daily, extremepre_daily, SHI, SPI,  type );
                
                % put all grids in a big matrix
                Com_hw_dr(i,j, type, :) = com_hw_dr(:,1) ;  Com_hw_ep(i,j,type, :) = com_hw_ep(:,1);
                Com_cw_dr(i,j,type,:) = com_cw_dr(:,1);  Com_cw_ep(i,j,type,:) = com_cw_ep(:,1);
                
%                 sum( com_hw_dr(:,1) )/72 % for check
%                 sum(~isnan(heatwave_daily) )/72
%                 sum(~isnan(drought_daily) )/72
%                 
            end
        end
    end
end
toc
% save compound_daily.mat  Com_hw_dr Com_hw_ep Com_cw_dr  Com_cw_ep 