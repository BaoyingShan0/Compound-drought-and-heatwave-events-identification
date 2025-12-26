function [SHI, h_ks, sim_cdf, emp_cdf, NSE_best, best_dist_idx, models, H_mk]=SHI_best_dist(Date, Data, scale, NSP)


% 2023.04.03

% advantages:
% 1. handle zero problem by conditional probability method
% 2. handle non stationarity problem by refereing to the past 30 years
% 3. handle fitting problem by offering several commonly used distributions for choosing the
% best according to the minimum AIC

% Input Data:
% Date: daily date time series, 3 columes should be [year, month, day]
% Data : daily Data vector not matrix, there 365 days for each year
% scale : days, like 3 days
% NSP: the variable for non stationarity. If NSP=30,  then consider the past 30 years as "normal" condition'
% if there is no input of NSP, that's mean data is stationary, then consider the whole date series as normal condition .

% output:
% SHI: daily SHI value % eg. if scale =3days, fist year data first 2 rows of SHI values are nan.
% h_ks: the KS test results for whether best distribution pass ks test at 0.05 significant level, if h=0, the distribution is OK.
% sim_cdf , emp_cdf: simulated and empirical cdf of Data
% NSE_best: Nash efficient value of the best distribution

% best_dist_idx: the index of best distribution for Data among models set in function
% models: all the distribution models for chosen
% H_mk: whether the XS is stationary or non stationariy, -1, 0, 1, using mk
% test at 0.05 significant level

% Na = nargin;
arguments %default value
    Date double = []
    Data double = []
    scale  double = 3
    NSP double = 30
end

%% ks.test.critical.value
% Load the critical value Dn0 and the corresponding significant level for the KS test: this table
% size is 2000*6 . the value of each columns is at
% significant 0.001,0.01,0.05,0.1,0.15,0.2, respectively
% rows mean the sample number from 1 to 5000, or larger
% calculated from the Dn0.R
% Dn0_level=[0.001,0.01,0.05,0.1,0.15,0.2];
sig_p=0.05;
load Dn0.mat

%% accumulative
calendar_days = 365; % There are 365 days for each calendar year.
years = Date(end,1) - Date(1,1) + 1; % Calculate the total number of years in this time series.
% accumulative climate variable
A1 = [];
for j = 1:scale
    A1 = [A1, Data(j:length(Data)-scale+j)];
end
XS = mean(A1, 2); % Calculate the mean of each row
XS = [nan(scale-1, 1); XS];

% Initialize variables
[SHI, h_ks, sim_cdf, emp_cdf, NSE_best, best_dist_idx] = deal(nan(length(Data), 1));
H_mk = nan(calendar_days,1);
alpha_mk=0.05;% alpha for mk ks test
%% SHI calculation

% for each calendar day
% for each yr
% if non stationary, past 30 years
% if stationary, no significant trend, then the expected condition is the long term mean; else, past 30 years, need to talk with bernard, niko


% tic
for j=1: calendar_days
    % j % indicate calendar day

    Xn=XS( j:calendar_days:length(XS) );
    [h,p, Z, h2] = mann_kendall(Xn,alpha_mk);
    H_mk(j)=h2;

    if h % if nonstationary

        i=1; % indicate year
        while i <= years
            % i

            if i==1
                tind=j: calendar_days: calendar_days*(NSP-1)+j; % for example, NSP=30 yeas, for the first 30 years in time series, normal condition is always the first 30 years,
                i=i+NSP-1;
            else % for the 31th years and after, the normal condition is the past 30 years, the (j, j-1, j-2, ..., j-29);
                tind=calendar_days*(i-30)+j : calendar_days: calendar_days*( i-1)+j;
            end

            Xn=XS(tind);
            if sum(Xn<0)==0 && sum(Xn==0)>0 % if all temperature is non negative and there are zeros, to avoild -inf in shi

                %             [SimCDF,ind_best, pd, rmsd_best, nse_best, models] = non_negtive_zero(Xn);
                SimCDF=nan( length(Xn), 1 ) ;
                [zeroa]=find(Xn==0);
                Xn_nozero=Xn; Xn_nozero(zeroa)=[];
                q=length(zeroa)/length(Xn);

                [ind_best, models] = best_dist(Xn_nozero);
                pd=fitdist( Xn_nozero , models(ind_best )  ) ; % parameters of the best distribution


                SimCDF(Xn~=0)=q+(1-q)*cdf(  pd, Xn_nozero ) ;
                SimCDF(Xn==0)=q;

                [resCDF,queryRes] = ecdf(Xn_nozero );
            else
                [ ind_best, models ]= best_dist(Xn) ;
                pd=fitdist( Xn , models(ind_best )  ) ; % parameters of the best distribution
                SimCDF=cdf(  pd, Xn );
                [resCDF,queryRes] = ecdf( Xn );
            end

            invcdf=icdf(pd, resCDF);
            queryRes=queryRes(2:end-1);  invcdf=invcdf(2:end-1); % no P=1 or P=0;
            delta=queryRes-invcdf;
            %  mean square deviation of variable
            %             rmsd_best= sqrt( sum(delta.^2)/ length(invcdf) ); %RMSE root
            nse_best=1- sum(delta.^2)/ sum(  ( invcdf-mean(invcdf) ) .^2) ;


            if i<=NSP
                SHI(tind, : )= norminv(SimCDF) ;
                sim_cdf(tind, : )=SimCDF;
                emp=nan(length(Xn),1);
                for ii=1: length(Xn)
                    emp(ii)= sum( Xn<=Xn(ii) )/length(Xn) ;
                end
                emp_cdf(tind, : )=emp;
                NSE_best( tind, : )=nse_best;
                %                 RMSD_best(tind, :)=rmsd_best;
                best_dist_idx( tind, : )=ind_best;
                % ks test
                h_ks(tind) =KsTest( sort(Xn), sig_p, pd, Dn0, Dn0_level);
            else
                SHI( tind(end), : )= norminv( SimCDF(end) ) ; %
                sim_cdf( tind(end), : )=SimCDF(end);
                emp_cdf(tind(end), : )=sum( Xn<=Xn(end) )/length(Xn);
                NSE_best( tind(end), : )=nse_best;
                best_dist_idx( tind(end), : )=ind_best;

                % ks test
                h_ks( tind(end) ) =KsTest(sort(Xn), sig_p, pd, Dn0, Dn0_level);
            end

            i=i+1;
        end

    else %else the data is stationary

        tind=j:calendar_days:length(XS);
        Xn=XS(tind);

        if sum(Xn<0)==0 && sum(Xn==0)>0 % if all temperature is non negative and there are zeros, to avoild -inf in shi

            %         [SimCDF,ind_best, pd, rmsd_best, nse_best, models] = non_negtive_zero(Xn);
            SimCDF=nan( length(Xn), 1 ) ;
            [zeroa]=find(Xn==0);
            Xn_nozero=Xn; Xn_nozero(zeroa)=[];
            q=length(zeroa)/length(Xn);

            [ind_best, models] = best_dist(Xn_nozero);
            pd=fitdist( Xn_nozero , models(ind_best )  ) ; % parameters of the best distribution


            SimCDF(Xn~=0)=q+(1-q)*cdf(  pd, Xn_nozero ) ;
            SimCDF(Xn==0)=q;
            [resCDF,queryRes] = ecdf(Xn_nozero );
        else
            [ ind_best, models ]= best_dist(Xn) ;
            pd=fitdist( Xn , models(ind_best )  ) ; % parameters of the best distribution
            SimCDF=cdf(  pd, Xn );
            [resCDF,queryRes] = ecdf( Xn );
        end

        SHI(tind, : )= norminv(SimCDF) ;
        sim_cdf(tind, : )=SimCDF;
        invcdf=icdf(pd, resCDF);
        queryRes=queryRes(2:end-1);  invcdf=invcdf(2:end-1); % no P=1 or P=0;
        delta=queryRes-invcdf;
        %  mean square deviation of variable
        %         rmsd_best= sqrt( sum(delta.^2)/ length(invcdf) ); %RMSE root
        nse_best=1- sum(delta.^2)/ sum(  ( invcdf-mean(invcdf) ) .^2) ;

        emp=nan(length(Xn),1);
        for ii=1: length(Xn)
            emp(ii)= sum( Xn<=Xn(ii) )/length(Xn) ;
        end
        emp_cdf(tind, : )=emp;
        NSE_best( tind, : )=nse_best;
        %         RMSD_best(tind, :)=rmsd_best;
        best_dist_idx( tind, : )=ind_best;
        % ks test
        h_ks(tind) =KsTest( sort(Xn), sig_p, pd, Dn0, Dn0_level);
    end

end
% toc
end
