function [SPI, h, sim_cdf , emp_cdf, NSE, best]=SPI_best_dist(Date, Data, scale, Model, NSP)
Na = nargin;
% 2021.10.22
% update the ks test method on 2022.1.22
% could handle zeros by conditional probability method
% pick the best distribution from several common used distributions
% according to AIC and ks test

% Input Data:
% Date: daily date time series, 3 columes should be [year, month, day]
% Data : daily Data vector not matrix, there 365 days for each year
% scale : days, like 15,30 days
% Model: given fitting distribution, for example, Model="gev". if Model="", then find the best distribution function among 8 models for data
% NSP: the variable for non stationarity. If NSP=30,  then consider the past 30 years as "normal" condition'
% if there is no input of NSP, that's mean data is stationary, then consider the whole date series as normal condition .

% output:
% daily SPI value % eg. if scale =15 days, fist year data first 14 rows of SHI values are nan.
% h: the KS test results for whether best distribution pass ks test at 0.05 significant level, if h=0, the distribution is OK.
% sim_cdf , emp_cdf: simulated and empirical cdf of Data
% NSE: Nash efficient value
% best: the index of best distribution for Data among model set in function
% best_fitting.m


% Load the critical value table: this table
% size is 2000*6 . the value of each columns is at
% significant 0.001,0.01,0.05,0.1,0.15,0.2, respectively
% rows mean the sample num from 1 to 2000
Dn0= readtable('Dn.csv', 'HeaderLines' ,1, 'ReadRowNames' , 1); % skip the first line: the variable name line.

%%  366 ---> 365 days
%For easier calculation of the drought index, all years are considered to have 365 days;
%for the leap years, 29 February is excluded and the rainfall on 28 February is recalculated as the average rainfall observed on 28 and 29 February.

mm = find( Date(:,2)==2 & Date(:,3)==29) ;
if ~isempty(mm)
    Data(mm-1,:)=(Data(mm-1,:)+Data(mm, :) )/2;
    Data(mm,:)=[];
end

% calculate SHI

%% “Best-fitting” picked based max AIC

days_year=365; % there are 365 days for each calendar year
years=Date(end,1)-Date(1,1)+1; %the total number of years in this time series.
% Data setting to scaled dataset
A1=[];
for is=1: scale,  A1=[A1, Data(is:length(Data)-scale+is)]; end
XS=sum(A1,2); % the sum
XS=[ nan(scale-1, 1); XS ];
[SPI,  h, sim_cdf, emp_cdf, NSE, best]=deal( nan( length(Data) , 1 ) );% deal

%%

% is=308; j=85; for test fitting

for is=1: days_year
     is
    if Na == 5 % if there is the input for nonstationary
        j=1;
        while j <= years
             j
            emp=[]; nse=nan; Best_fitting_H=nan;
            
            if j==1
                tind=is: days_year: days_year*(NSP-1)+is; % for example, NSP=30 yeas, for the first 30 years in time series, normal condition is always the first 30 years,
                j=j+NSP-1;
            else % for the 31th years and after, the normal condition is the past 30 years, the (j, j-1, j-2, ..., j-29);
                tind=days_year*(j-30)+is : days_year: days_year*( j-1)+is;
            end
            
            Xn=XS(tind); SimCDF=nan( length(Xn), 1 ) ; 
            %             best_fitting(Xn, 'p')
            % replace the nan with average value, adding by Baoying 2021.8.12
            % [~, bb]=rmmissing(Xn);
            % Xn(bb)=nanmean(Xn); % use mean to replace the missing value, this does not commonly happen
            [zeroa]=find(Xn==0);
            Xn_nozero=Xn; Xn_nozero(zeroa)=[];
            q=length(zeroa)/length(Xn);
            
            if Model=="" % 如果 model 输入为空
                %  find the "best-fitting" distribution and this best dist
                %  need to pass ks test
                 [ Best_fitting_H, nse, model ]= best_fitting( Xn_nozero, '', true, "AIC") ; %AIC and kstest as criteria
               best_fitting(data, for_plot, for_ks, for_criteria)
                 pd=fitdist( Xn_nozero , model(Best_fitting_H )  ) ; % parameters of the best distribution
            else
                pd=fitdist( Xn_nozero, Model ) ;% parameters of the given distribution
                nse = Nash( Xn_nozero, pd);
            end
            
            %             SimCDF=cdf(  pd, Xn );
            SimCDF(Xn~=0)=q+(1-q)*cdf(  pd, Xn_nozero ) ;
            SimCDF(Xn==0)=q;
            
            if j<=NSP
                SPI(tind, : )= norminv(SimCDF) ; %SHI
                sim_cdf(tind, : )=SimCDF;
                for i=1: length(Xn)
                    emp(i)= sum( Xn<=Xn(i) )/length(Xn) ;
                end
                emp_cdf(tind, : )=emp;
                NSE( tind, : )=nse;
                
                % ks test
                sig_p=0.05;
                h(tind) =ksTestby( sort(Xn), sig_p, pd, Dn0);
            else
                SPI( tind(end), : )= norminv( SimCDF(end) ) ; %
                sim_cdf( tind(end), : )=SimCDF(end);
                emp_cdf(tind(end), : )=sum( Xn<=Xn(end) )/length(Xn);
                NSE( tind(end), : )=nse;
                best( tind(end), : )=Best_fitting_H;
                % ks test
                sig_p=0.05;
                h( tind(end) ) =ksTestby(sort(Xn), sig_p, pd, Dn0);
                
            end
            j=j+1;
        end
        
    else % the data is stationary
        
        tind=is:days_year:length(XS);
        Xn=XS(tind);  SimCDF=nan( length(Xn), 1 );
        % is=4
        %best_fitting(Xn(14:43), 'p')
        % replace the nan with average value, adding by Baoying 2021.8.12
        % [~, bb]=rmmissing(Xn);
        % Xn(bb)=nanmean(Xn); % use mean to replace the missing value, this does not commonly happen
        
        %% how to handle 0? 2022.1.21
        [zeroa]=find(Xn==0);
        Xn_nozero=Xn;Xn_nozero(zeroa)=[];
        q=length(zeroa)/length(Xn);
        %     parm=gamfit(Xn_nozero);
        %     Gam_xs=q+(1-q)*gamcdf(Xn,parm(1),parm(2));
        %     Z(tind)=norminv(Gam_xs);
        if Model ==""  % find the "best-fitting" distribution && this best pass ks test
            [ Best_fitting_H, nse, model ]= best_fitting( Xn_nozero, '', true, "AIC") ; %AIC
            pd=fitdist( Xn_nozero , model(Best_fitting_H )  ) ; % parameters of the best distribution
            best( tind, : )=Best_fitting_H;
        else
            pd=fitdist(Xn_nozero, Model) ;% parameters of the given Model distribution
            nse=Nash(Xn_nozero, pd);
        end
        SimCDF(Xn~=0)=q+(1-q)*cdf(  pd, Xn_nozero );
        SimCDF(Xn==0)=q;
        %         if Model ==""  % find the "best-fitting" distribution && this best pass ks test
        %             [ Best_fitting_H, nse, model ]= best_fitting( Xn, '', true ) ;
        %
        %             pd=fitdist( Xn , model(Best_fitting_H )  ) ; % parameters of the best distribution
        %             best( tind, : )=Best_fitting_H;
        %         else
        %             pd=fitdist(Xn, Model) ;% parameters of the given Model distribution
        %             nse=Nash(Xn, pd);
        %         end
        %         SimCDF=cdf(  pd, Xn );
        
        SPI(tind, : )=norminv(SimCDF);
        sim_cdf( tind, : )=SimCDF;
        NSE( tind, : )=nse;
        
        for i=1: length(Xn)
            emp(i)=sum( Xn<=Xn(i) ) / length(Xn);
        end
        emp_cdf( tind, : )=emp;
        
        sig_p=0.05;
        h(tind) = ksTestby(sort(Xn), sig_p, pd, Dn0) ; % could ignore ks test
    end
    
end


%%  子函数
function [h, testResults] =ksTestby(data, p, pd ,Dn0)

% 2022.1.11
% this function is to test if data follow give model diatribution by ks test
% at p significant level

% data: one colunm is one series of data
% p: significant level
% model: given distribution model, like "gev", "gpd"

N =size(data,2);

ps=[0.001	0.01	0.05	0.1	0.15	0.2];
Dn=Dn0(:, ps==p);

for i =1:N
    % exclude nan value
    x = data(:,i);
    x(isnan(x))=[];
    
    % get theritical cdf based on given pd
    xx=sort(x) ;
    theop = cdf(pd, xx) ;
    
    % empirical
    n=length(x);
    obspu = [1:1: n]/n;
    obspl =[0:1:n-1] / n;
    ks = max( max( abs(theop - obspu'), abs(theop - obspl') ) );
    
    dn=Dn{n,1};
    if ks <= dn
        testResults{i}="not rejected"; h(i)=0;
    else
        testResults{i}="rejected"; h(i)=1;
    end
    
end

