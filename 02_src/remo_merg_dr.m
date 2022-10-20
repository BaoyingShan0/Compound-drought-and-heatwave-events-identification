function [best_combination_for_drought, Table_dr_RM, criteria]=remo_merg_dr(Date, SPI, REMO_d, MERG_d, p, start_th_d, end_th_d)

% 2022.2.12,
% changed on 2022.7.13

%
% Aim: this function could get cases where their severity follows GEV and
% inter-arrival time follow exponential distribution (by ks test)
%  and given the best one which has maximum events number.
% However, it may happen that the two criteria are too strict to see permutations.
% In that case, only AIC is used as a criterion.

% this function could get the final  thresholds for removing and
% merging and don't need to go the R and excel. Why not use R? the speed
% for fitting in matlab is quicker and now don't need transfer dataset between
% two software.

% The start and end thresholds for pre identifying is -1 for drought and 1 for heatwave in this fuction.

% input:
% Date: [year, month, day]; 365 days each year
% SPI, SHI: daily value, correpondeing with date
% REMO_d, MERG_d,  REMO_h, MERG_h: possible threholds for removing and
% merging, eg REMO_d=[45:-1:0]; MERG_d=[15:-1:0]; REMO_h=[15:-1:0]; MERG_h=[15:-1:0];
% p: significant level for ks test
% start_th_d, end_th_d, start_th_h, end_th_h: start and thresholds for dr and hw

% output:
% best_combinations: 2*2 matrix,  [remove threshold for dr,  merge threshold for dr; remove threshold for hw,  merge threshold for hw]
% Table_dr_RM and Table_hw_RM: the table including all results: cases, best distribution for severity,
% ks test if arrival follow exp, annual days, annual events.

% based on remo_merge, seperate the drought and heatwave
% nargin

if nargin < 5
    p=0.05;
end
if nargin<7
    start_th_d=-1; end_th_d=start_th_d;
end

%% I  AIC + RMSD
criteria=["AIC & RMSD"];

years=size(Date,1)/365;
who_first='REMO';

CASES_dr = [];
for a=1: length(REMO_d)
    for b=1: length(MERG_d)
        CASES_dr=[CASES_dr; REMO_d(a),  MERG_d(b)];% removing, merging of drought,
    end
end

N= size(CASES_dr,1);
[annual_events_dr, annual_days_dr, ind_dr, h_dr, best_dist_dr]=deal( nan( N , 1 ) );% deal
drought_daily=[]; best_dist_dr=string( best_dist_dr );
for j=1: N
%     j
    if who_first=='REMO'
        drought_daily=PRM2_drought_identification( Date, SPI, start_th_d,end_th_d ,...
            CASES_dr(j,1),  CASES_dr(j,2) );  %xlim( time_lim ); datetick('x','yyyy-mm' );
        % drought_daily=PRM2_drought_identification( Date, SPI, start_th_d,end_th_d, 19, 0);
    elseif who_first=='MERG'
        % try  remove first then remove
        drought_daily=PMR2_drought_identification( Date, SPI(:,k), start_th_d,end_th_d ,...
            CASES_dr(j,1),  CASES_dr(j,2) );
    end
    size_dr=size(drought_daily,2);
    M= max( drought_daily(:,size_dr) ) ;
    dr= nan(M, 5);
    if M>10   %if the drought numbers is less than 10, then skip
        for i=1: M%
            aa=drought_daily( drought_daily(:, size_dr )==- i,:) ;
            dr(i,1:4)=[i, size(aa,1), -sum(aa(:,4) ), -sum(aa(:,4) )/size(aa,1)] ; %num, duration, severity, intensity
            if i<2
                dr(i,5)=nan;
            else
                bb=find(drought_daily(:, size_dr )==-i  ) ;
                cc=find( drought_daily(:, size_dr )==-i+1 ); % last events
                dr(i,5)=bb(1)-cc(end)-1; % the duration of the non-dry period i
            end
        end
        
        arrivals_dr=dr(2:end, 2)+dr(2:end, 5);
        
        annual_days_dr ( j, : )=sum( dr ( :,2) )/ years; % annual days in drought
        annual_events_dr ( j, : )=size(dr,1)/years; % annual drought events
        
        % is the best distribution for severity GEV? based on AIC
        [ind_dr(j), ~, models]=best_fitting( dr(:, 3), "", true );
        best_dist_dr(j)=  models( ind_dr(j) );
        
        % kstest to get if inter arrival time follow exponential
        h_dr(j) =ksTestby(arrivals_dr, p, "exp"); % if h_dr==0,  arrivals follow gev
        
    end
end

% the index where severity follows gev and arrivals follow exponential distribution
AA= best_dist_dr=="gev" & h_dr==0 & annual_days_dr<365/2 ;
% the index where severity follows gev and arrivals follow exponential
% distribution ( h_dr==0 ); and  annual days in drought should less than half a year, a very wide constraints
% aa= best_dist_dr=="gev" & h_dr ==0 ;

%% II: if there is no combination left using min AIC + min RMSD, then only use AIC

if sum(AA)==0 %
    criteria=["AIC"];
    N= size(CASES_dr,1);
    [annual_events_dr, annual_days_dr, ind_dr, h_dr, best_dist_dr]=deal( nan( N , 1 ) );% deal
    drought_daily=[]; best_dist_dr=string( best_dist_dr );
    for j=1: N
%         j
        if who_first=='REMO'
            drought_daily=PRM2_drought_identification( Date, SPI, start_th_d,end_th_d ,...
                CASES_dr(j,1),  CASES_dr(j,2) );  %xlim( time_lim ); datetick('x','yyyy-mm' );
            % drought_daily=PRM2_drought_identification( Date, SPI, start_th_d,end_th_d, 19, 0);
        elseif who_first=='MERG'
            % try  remove first then remove
            drought_daily=PMR2_drought_identification( Date, SPI(:,k), start_th_d,end_th_d ,...
                CASES_dr(j,1),  CASES_dr(j,2) );
        end
        size_dr=size(drought_daily,2);
        M= max( drought_daily(:,size_dr) ) ;
        dr= nan(M, 5);
        if M>10   %if the drought numbers is less than 10, then skip
            for i=1: M%
                aa=drought_daily( drought_daily(:, size_dr )==- i,:) ;
                dr(i,1:4)=[i, size(aa,1), -sum(aa(:,4) ), -sum(aa(:,4) )/size(aa,1)] ; %num, duration, severity, intensity
                if i<2
                    dr(i,5)=nan;
                else
                    bb=find(drought_daily(:, size_dr )==-i  ) ;
                    cc=find( drought_daily(:, size_dr )==-i+1 ); % last events
                    dr(i,5)=bb(1)-cc(end)-1; % the duration of the non-dry period i
                end
            end
            
            arrivals_dr=dr(2:end, 2)+dr(2:end, 5);
            
            annual_days_dr ( j, : )=sum( dr ( :,2) )/ years; % annual days in drought
            annual_events_dr ( j, : )=size(dr,1)/years; % annual drought events
            
            % is the best distribution for severity GEV? only based on AIC
            [ind_dr(j), ~, models]=best_fitting(  dr(:, 3), "", true, "AIC"  );
            best_dist_dr(j)=  models( ind_dr(j) );
            
            % kstest to get if inter arrival time follow exponential
            h_dr(j) =ksTestby(arrivals_dr, p, "exp"); % if h_dr==0,  arrivals follow gev
            
        end
    end
    
    AA= best_dist_dr=="gev" & h_dr==0 & annual_days_dr<365/2 ;
end

%% final output
% table results for output
reject_or_not= [ "not rejected"; "rejected"];
testResults= string( nan( N , 1 ) ); cc=~isnan(h_dr);
testResults(cc) = reject_or_not( h_dr(cc)+1) ;
Table_dr_RM = table(CASES_dr(:,1), CASES_dr(:,2), best_dist_dr, testResults,  annual_days_dr, annual_events_dr );
Table_dr_RM.Properties.VariableNames(1: 6)= {'Threshold_REMO', 'Threshold_MERG', ...
    'Best distribution for severity', 'ks test if arriavl follows exponential', 'Annual days in drought' , 'annual of events'};

annual_events_dr(~AA)=nan;
[~, final_ind_dr]=nanmax(annual_events_dr);
best_combination_for_drought=table2array( Table_dr_RM(final_ind_dr, :) );