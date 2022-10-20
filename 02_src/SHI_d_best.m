function [SHI, h, sim_cdf , emp_cdf, NSE, best]=SHI_d_best(Date, Data, scale, Model, NSP)

% 2021.10.22
% update the ks test method on 2022.1.11
% Standardized Heatwave Index, similar with SPI
% pick the best distribution from several common used distributions

% Input Data:
% Date: daily date time series, 3 columes should be [year, month, day]
% Data : daily Data vector not matrix, there 365 days for each year
% scale : days, like 3 days, or 5 days or 7 days
% Model: given fitting distribution, for example, Model="gev". if Model="", then find the best distribution function among 8 models for data
% NSP: is the period value for nonstationary, if the Data is nonstationary
% then consider the past 30 years as "normal" condition, instead of the
% whole Date series.

% output:
% SHI daily SHI value % eg. if scale =3 days, fist year data first 2 rows of SHI values are nan.
% h: the KS test value, if h=0, the distribution is OK.


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
    Data(mm-1,4: 6)=(Data(mm-1,:)+Data(mm, :) )/2;
    Data(mm,:)=[];
end

% calculate SHI

%% “Best-fitting” picked based max AIC

days_year=365; % there are 365 days for each calendar year
years=Date(end,1)-Date(1,1)+1; %the total number of years in this time series.
% Data setting to scaled dataset
A1=[];
for is=1: scale,  A1=[A1, Data(is:length(Data)-scale+is)]; end
XS=mean(A1,2); % the mean DMT over past scale days, (T + Ti-1 + Ti-2)/3
XS=[ nan(scale-1, 1); XS ];
[SHI,  h, sim_cdf, emp_cdf, NSE, best]=deal( nan( length(Data) , 1 ) );% deal

%%

% is=308; j=85; for test fitting

for is=1: days_year
    is
    if nargin == 5 % if there is the input for nonstationary
        j=1;
        while j <= years
            j;
            emp=[]; nse=nan; Best_fitting_H=nan;
            
            if j==1
                tind=is: days_year: days_year*(NSP-1)+is; % for example, NSP=30 yeas, for the first 30 years in time series, normal condition is always the first 30 years,
                j=j+NSP-1;
            else % for the 31th years and after, the normal condition is the past 30 years, the (j, j-1, j-2, ..., j-29);
                tind=days_year*(j-30)+is : days_year: days_year*( j-1)+is;
            end
            
            Xn=XS(tind);
            
            if sum(Xn<0)==0 && sum(Xn==0)>0 % if all temperature is non negative and tere are zeros, to avoild -inf in shi
                SimCDF =non_negtive_zero( Xn, Model );
            else
                % replace the nan with average value, adding by Baoying 2021.8.12
                % [~, bb]=rmmissing(Xn);
                % Xn(bb)=nanmean(Xn); % use mean to replace the missing value, this does not commonly happen
                
                if Model=="" % 如果 model 输入为空
                    %  find the "best-fitting" distribution and this best dist
                    %  need to pass ks test
                    [ Best_fitting_H, nse, model ]= best_fitting(Xn, '', true, "AIC") ;
                    pd=fitdist( Xn , model(Best_fitting_H )  ) ; % parameters of the best distribution
                else
                    pd=fitdist(Xn, Model ) ;% parameters of the given distribution
                    nse = Nash( Xn, pd);
                end
                
                SimCDF=cdf(  pd, Xn );
            end
            
            if j<=NSP
                SHI(tind, : )= norminv(SimCDF) ; %SHI
                sim_cdf(tind, : )=SimCDF;
                for i=1: length(Xn)
                    emp(i)= sum( Xn<=Xn(i) )/length(Xn) ;
                end
                emp_cdf(tind, : )=emp;
                NSE( tind, : )=nse;
                
                % ks test
                sig_p=0.05;
                h(tind) =ksTestby(sort(Xn), sig_p, pd, Dn0);
            else
                SHI( tind(end), : )= norminv( SimCDF(end) ) ; %SHI
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
        Xn=XS(tind);
        
        if sum(Xn<0)==0 && sum(Xn==0)>0 % if all temperature is non negative and tere are zeros
            SimCDF =non_negtive_zero( Xn, Model );
        else
           % [~, bb]=rmmissing(Xn);
            % Xn(bb)=nanmean(Xn); % use mean to replace the missing value, this does not commonly happen
            
            if Model ==""  % find the "best-fitting" distribution && this best pass ks test
                [ Best_fitting_H, nse, model ]= best_fitting( Xn, '', true, "AIC" ) ; %AIC  as criteria
                pd=fitdist( Xn , model(Best_fitting_H )  ) ; % parameters of the best distribution
                best( tind, : )=Best_fitting_H;
            else
                pd=fitdist(Xn, Model) ;% parameters of the given Model distribution
                nse=Nash(Xn, pd);
            end
            
            SimCDF=cdf(  pd, Xn );
        end
        SHI(tind, : )=norminv(SimCDF);
        sim_cdf( tind, : )=SimCDF;
        NSE( tind, : )=nse;
        
        for i=1: length(Xn)
            emp(i)=sum( Xn<=Xn(i) ) / length(Xn);
        end
        emp_cdf(tind, : )=emp;
        
        sig_p=0.05;
        h(tind) = ksTestby(sort(Xn), sig_p, pd, Dn0) ;
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
end



function [SimCDF] =non_negtive_zero(Xn, Model)
        % if all temperature is non negative and there are zeros
        
        SimCDF=nan( length(Xn), 1 ) ;
        % replace the nan with average value, adding by Baoying 2021.8.12
        %[~, bb]=rmmissing(Xn);
        %Xn(bb)=nanmean(Xn); % use mean to replace the missing value, this does not commonly happen
        [zeroa]=find(Xn==0);
        Xn_nozero=Xn; Xn_nozero(zeroa)=[];
        q=length(zeroa)/length(Xn);
        
        if Model=="" % 如果 model 输入为空
            %  find the "best-fitting" distribution and this best dist
            %  need to pass ks test
            [ Best_fitting_H, nse, model ]= best_fitting(Xn_nozero, '', true, "AIC" ) ;
            pd=fitdist( Xn_nozero , model(Best_fitting_H )  ) ; % parameters of the best distribution
        else
            pd=fitdist( Xn_nozero, Model ) ;% parameters of the given distribution
            nse = Nash( Xn_nozero, pd);
        end
        
        SimCDF(Xn~=0)=q+(1-q)*cdf(  pd, Xn_nozero ) ;
        SimCDF(Xn==0)=q;
end

end

%
% % 子函数
% function [h, testResults] =ksTestby1(data, p, model)
%
% % 2022.1.11
% % this function is to test if data follow give model diatribution by ks test
% % at p significant level
%
% % data: one colunm is one series of data
% % p: significant level
% % model: given distribution model, like "gev", "gpd"
%
% N =size(data,2);
%
% ps=[0.001	0.01	0.05	0.1	0.15	0.2];
% Dn=Dn(:, ps==p);
%
% for i =1:N
%     % exclude nan value
%     x = data(:,i);
%     x(isnan(x))=[];
%
%
%     % fit distribution and get theritical cdf
%     pd=fitdist(x, model ) ;
%     xx=sort(x) ;
%     theop = cdf(pd, xx) ;
%
%     % empirical
%     n=length(x);
%     obspu = [1:1: n]/n;
%     obspl =[0:1:n-1] / n;
%     ks = max( max( abs(theop - obspu'), abs(theop - obspl') ) );
%
%     dn=Dn{n,1};
%     if ks <= dn
%         testResults{i}="not rejected"; h(i)=0;
%     else
%         testResults{i}="rejected"; h(i)=1;
%     end
%
% end
%
%
% ks
%


%{

% %%
% days_year=365;  years=Date(end,1)- Date(1,1)+1;
% reshaped_DMT=reshape( DMT, 365, [])';
% DDMT=reshaped_DMT;
% %  Best_fitting_T= nan(years, days_year);  NSE_T= nan(years, days_year);
% NSE_T_gev=nan(years, days_year);
% SHI= nan(years, days_year);
% tic
% % double loop makes the running slower:  20 minutes
% for i=1: days_year
%
%     for j=1: years
%
%         if j<=30
%             data=DDMT(1:30 , i);
%         else
%             data=DDMT( j-29 : j , i); % talk with niko, this is ok
%         end
%
% %         [ Best_fitting_T( j, i ), NSE_T( j, i), model ]= best_fitting( data ) ;
%          pd=fitdist( data, "gev"  ) ; % parameters of the given distribution
%
%          SimCDF=cdf(  pd, DDMT(j , i) ); % "relative heatwave": the relative is compared to the "normal", and the "normal" is defined as the past 30 years
%          SHI( j, i )=norminv(SimCDF); %SHI
%
%          [NSE_T_gev(j,i), ~]= nse_distribution( data , [ "gev"]) ; %NSE
%
%          % KS test
%          [resCDF, queryRes] = ecdf( data );
%          pd= fitdist( data, "gev") ;
%          invcdf=icdf(pd, resCDF(2:end-1) );
%          [ h( j, i ), p( j, i ) ] = kstest2( queryRes(2:end-1),  invcdf ); %h=0, gev is ok
%
%     end
%
% end
% toc
% sum(sum(h))
% NSE_T_gev0=NSE_T_gev; SHI0=SHI;
% NSE_T_gev=reshape(NSE_T_gev', [], 1); % Best_fitting_T= reshape ( Best_fitting_T',[],1 ) ;
% SHI= reshape(SHI0',[],1 );
%
% %%
% % don't consider the influence of nonstationary
%
% Best_fitting_T_o= nan(1, days_year);  NSE_T_o= nan(1, days_year);
% SHI_o= nan(years, days_year);
% for i=1: days_year
%        % [ Best_fitting_T_o(1, i ), NSE_T_o( 1, i), model ]= best_fitting( DDMT(:,i) ) ;
%         % pd=fitdist( DDMT( : , i ) , model( Best_fitting_T_o( 1, i ) )  ) ; % parameters of the best distribution
%         pd=fitdist( DDMT( : , i ) , "gev"  ) ;
%         SimCDF=cdf(  pd, DDMT(: , i)  ); % "relative heatwave": the relative is compared to the "normal", and the "normal" is defined as the past 30 years
%         SHI_o(:, i )=norminv(SimCDF);
% end
% SHI_o= reshape(SHI_o',[],1 );

%}
