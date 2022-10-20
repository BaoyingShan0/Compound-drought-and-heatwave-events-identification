function [best_combination_for_heatwave, Table_hw_RM, criteria]=remo_merg_hw(Date, SHI,   REMO_h, MERG_h, p, start_th_h, end_th_h)

% similar with remo_merg_dr
% 2022.7.1: 

if nargin < 5
    p=0.05;
end
if nargin<7
    start_th_h=1; end_th_h=start_th_h;
end

years=size(Date,1)/365;
who_first='REMO';

CASES_hw=[];
for a=1:length(REMO_h)
    for b=1:length(MERG_h)
        CASES_hw=[CASES_hw; REMO_h(a),  MERG_h(b)];% removing, merging of drought,
    end
end

%% I  AIC + RMSD
criteria=["AIC & RMSD"];

N= size(CASES_hw,1);
[annual_events_hw, annual_days_hw, ind_hw, h_hw, best_dist_hw]=deal( nan( N , 1 ) );% deal
heatwave_daily=[]; best_dist_hw=string( best_dist_hw );
for j=1: N
%     j
    % removing first
    heatwave_daily=PRM2_heatwave_identification( Date, SHI, start_th_h, end_th_h , CASES_hw(j,1),  CASES_hw(j,2) );
    size_hw=size(heatwave_daily,2);
    
    M = max( heatwave_daily(:,size_hw) ) ;
    hw = nan( M, 5);
    if M>10   %if the drought numbers is less than 10, then skip
        for i=1: M%
            aa=heatwave_daily( heatwave_daily(:, size_hw )==-i,:) ;
            hw(i,1:4)=[i, size(aa,1),  sum(aa(:,4) ), sum(aa(:,4) )/size(aa,1)] ; %num, duration, severity, intensity
            if i<2
                hw(i,5)=nan;
            else
                bb=find(heatwave_daily(:, size_hw )==-i  ) ;
                cc=find( heatwave_daily(:, size_hw )==-i+1 ); % last events
                hw(i,5)=bb(1)-cc(end)-1; % the duration of the non-dry period i
            end
        end
        
        arrivals_hw=hw(2:end, 2)+hw(2:end, 5);
        annual_events_hw ( j, : )=size(hw,1)/years; % annual drought events
        annual_days_hw(j, :)=sum(hw(:,2) )/ years; %annual days in drought
        
        % is the best distribution for severity GEV? based on AIC +RMSD
        [ ind_hw(j), ~, models]=best_fitting( hw(:, 3), "", true  );
        best_dist_hw(j) = models( ind_hw(j) );
        
        % kstest to get if inter arrival time follow exponential
        h_hw(j) =ksTestby(arrivals_hw, p, "exp"); % if h_dr==0,  arrivals follow gev
    end
end
% the index where severity follows gev and arrivals follow exponential distribution
AA= best_dist_hw=="gev" & h_hw==0 & annual_days_hw<365/2 ;

%% II: if there is no combination left, only AIC
if sum(AA)==0 %
    criteria=["AIC"];
    N= size(CASES_hw,1);
    [annual_events_hw, annual_days_hw, ind_hw, h_hw, best_dist_hw]=deal( nan( N , 1 ) );% deal
    heatwave_daily=[]; best_dist_hw=string( best_dist_hw );
    for j=1: N
%         j
        % removing first
        heatwave_daily=PRM2_heatwave_identification( Date, SHI, start_th_h, end_th_h , CASES_hw(j,1),  CASES_hw(j,2) );
        size_hw=size(heatwave_daily,2);
        
        M = max( heatwave_daily(:,size_hw) ) ;
        hw = nan( M, 5);
        if M>10   %if the drought numbers is less than 10, then skip
            for i=1: M%
                aa=heatwave_daily( heatwave_daily(:, size_hw )==-i,:) ;
                hw(i,1:4)=[i, size(aa,1),  sum(aa(:,4) ), sum(aa(:,4) )/size(aa,1)] ; %num, duration, severity, intensity
                if i<2
                    hw(i,5)=nan;
                else
                    bb=find(heatwave_daily(:, size_hw )==-i  ) ;
                    cc=find( heatwave_daily(:, size_hw )==-i+1 ); % last events
                    hw(i,5)=bb(1)-cc(end)-1; % the duration of the non-dry period i
                end
            end
            
            arrivals_hw=hw(2:end, 2)+hw(2:end, 5);
            annual_events_hw ( j, : )=size(hw,1)/years; % annual drought events
            annual_days_hw(j, :)=sum(hw(:,2) )/ years; %annual days in drought
            
            % is the best distribution for severity GEV? only based on AIC
            [ ind_hw(j), ~, models]=best_fitting( hw(:, 3), "", true, "AIC"  );
            best_dist_hw(j) = models( ind_hw(j) );
            
            % kstest to get if inter arrival time follow exponential
            h_hw(j) =ksTestby(arrivals_hw, p, "exp"); % if h_dr==0,  arrivals follow gev
        end
    end
    AA= best_dist_hw=="gev" & h_hw==0 & annual_days_hw<365/2 ;
end

%% final output

% table results for output
reject_or_not = [ "not rejected"; "rejected"];
testResults_hw = string( nan( N , 1 ) ); cc = ~isnan(h_hw);
testResults_hw(cc) = reject_or_not( h_hw(cc)+1) ;

Table_hw_RM = table(CASES_hw(:,1), CASES_hw(:,2), best_dist_hw,  testResults_hw,  annual_days_hw, annual_events_hw );
Table_hw_RM.Properties.VariableNames(1: 6)= {'Threshold_REMO', 'Threshold_MERG', ...
    'Best distribution for severity', 'ks test if arriavl follows exponential', 'Annual days in heatwave' , 'annual of events'};

% get the final combination based on maximum annual events number
annual_events_hw(~AA)=nan;
[~, final_ind_hw]=nanmax(annual_events_hw);
best_combination_for_heatwave=table2array( Table_hw_RM(final_ind_hw, :) );

