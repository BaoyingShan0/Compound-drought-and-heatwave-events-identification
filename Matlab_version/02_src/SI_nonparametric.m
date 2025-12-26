function [SI]=SI_nonparametric(Date, Data, scale, NSP)

% use nonparametric method to calculate the standardized index
% from: Farahmand, Alireza, and Amir AghaKouchak. "A Generalized Framework for Deriving Nonparametric Standardized Drought Indicators." Advances in Water Resources 76 (February 1, 2015): 140–45.
% https://doi.org/10.1016/j.advwatres.2014.11.012.
% 2024.11.11, by Baoying Shan

% Input Data:
% Date: daily date time series, 3 columes should be [year, month, day]
% Data : daily Data vector not matrix, there 365 days for each year
% scale : days, like 15,30 days
% NSP: the variable for non stationarity. If NSP=30,  then consider the past 30 years as "normal" condition'
% if there is no input of NSP, that's mean data is stationary, then consider the whole date series as normal condition .

% output:
% SI: daily standardized index value % eg. if scale =15 days, fist year data first 14 rows of SI values are nan.


% Na = nargin;
arguments %default value
    Date double = []
    Data double = []
    scale  double = 15
    NSP double = 30
end


%% accumulative
calendar_days = 365; % There are 365 days for each calendar year.
years = Date(end,1) - Date(1,1) + 1; % Calculate the total number of years in this time series.
% accumulative climate variable
A1 = [];
for j = 1:scale
    A1 = [A1, Data(j:length(Data)-scale+j)];
end
XS = sum(A1, 2); % Calculate the sum of each row
XS = [nan(scale-1, 1); XS];

% Initialize variables
SI = nan(length(Data), 1);% RMSD_best,
H_mk = nan(calendar_days,1);
alpha_mk=0.05;% alpha for mk ks test

%% SI calculation
% for each calendar day
% for each yr
% if non stationary, past 30 years
% if stationary, no significant trend, then the expected condition is the long term mean; else, past 30 years, need to talk with bernard, niko

for j=1: calendar_days
    % day = j % indicate calendar day
    Xn=XS(j:calendar_days:length(XS) );
    [h,~, ~, h2] = mann_kendall(Xn,alpha_mk);
    H_mk(j)=h2;

    if h % if there is the input for nonstationary
        i=1; % indicate year
        while i <= years
            if i==1
                tind=j: calendar_days: calendar_days*(NSP-1)+j; % for example, NSP=30 yeas, for the first 30 years in time series, normal condition is always the first 30 years,
                i=i+NSP-1;
            else % for the 31th years and after, the normal condition is the past 30 years, the (j, j-1, j-2, ..., j-29);
                tind=calendar_days*(i-30)+j : calendar_days: calendar_days*( i-1)+j;
            end

            Xn=XS(tind);
            % same as the SDAT by Farahmand in 2015
            n=length(Xn);
            ii=zeros(n,1);
            for tt=1:n
                ii(tt,1)=sum(Xn(:,1)<=Xn(tt,1)); % the nuber of elements in the sample <= Xn(tt)
            end
            q=(ii-0.44)./(n+0.12);

            if i<=NSP
                SI(tind, : ) = norminv(q) ;
            else
                SI(tind(end), : )= norminv( q(end) ) ; %
            end

            i=i+1;
        end
    else %else the data is stationary
        tind=j:calendar_days:length(XS);
        Xn=XS(tind);
        n=length(Xn);
        ii=zeros(n,1);
        for tt=1:n
            ii(tt,1)=sum(Xn(:,1)<=Xn(tt,1)); % the nuber of elements in the sample <= Xn(tt)
        end
        q=(ii-0.44)./(n+0.12);
        SI(tind)=norminv(q);
    end

end

