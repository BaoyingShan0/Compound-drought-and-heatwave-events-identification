
function [identification_results]=PRM_extreme_identification(extreme_type, Date, extreme_index0, start_th, end_th, REMO, MERG)

%  2024.01.18
% Outputs: identification_results: each column is [ year, month, day, drought index, dry or non dry
% symbol after pre identification (1 means in dry period and 0 is in non
% dry), period orders after preidentification (negative is the order of
% dry period and positive is for non-dry period), dry or non dry symbol
% after removing, period orders after removing, dry or non dry symbol
% after merging, period orders after merging]

% Inputs:
%Date: time: 3 columns [year, month, day];
% extreme_index0: standardized extreme index value, like SPI or SHI, time series;
% STEP 1: pre identification to get dry period and non dry period. a dry period starts when index is less or equal than the start_th
% and end when the index becomes larger than end_th;
%start_th and end_th can be the same value.

% STEP 2: remove or delete minor hot period
% if the duration of the dry period is samller than REMO_time, this dry period is not considered as drought.

% STEP 3: merge the non-independent near hot periods
% if the sum of (start_th (1)-SHI ) during intevals is less then  MERG,
% then merge the two near hot-spells two one

%% 0. Extreme type
if ismember(extreme_type, ["dr","cw","d","c"]) 
    index = - extreme_index0;
    start_th = -start_th;
    end_th = -end_th;
elseif ismember(extreme_type, ["wet","hw","p","h"])
    index = extreme_index0;
end

%% 1. pre-identification
N=size(index,1 );
identification_results=[Date, index, nan(N,1) ];

for n=1:N
    if n==1
        if identification_results(n,4)>=start_th,  identification_results( n, 5 )=1;
        else, identification_results( n, 5 )=0; end
    else
        if identification_results(n,4)>=start_th
            identification_results( n, 5 )=1; % 1 means in dry period
        elseif identification_results( n-1, 5 )==1 &&  identification_results(n,4) > end_th
            identification_results( n, 5 )=1; % 1 means in dry period
        else
            identification_results( n, 5 )=0; %0 means in non dry period
        end
    end
end

identification_results( :,6 )=nan(N,1);
aa=0; % the order of non dry period
bb=0; % the order of dry period

% the number of dry periods and wet periods
% all the days in the same period have same number
for n=1:N
    if n==1
        if identification_results( n, 5 )==1
            identification_results( n, 6 )=bb;
        elseif identification_results( n, 5 )==0
            identification_results( n, 6 )=aa;
        end
    else
        if identification_results( n, 5 )==1 && identification_results( n-1, 5 )==1
            identification_results( n, 6 )=bb;
        elseif identification_results( n, 5 )==1 && identification_results( n-1, 5 )==0
            bb=bb-1;
            identification_results( n, 6 )=bb;  % the order of dry period, bb=-2 means the second dry period from the begaining
        elseif identification_results( n, 5 )==0 && identification_results( n-1, 5 )==1
            aa=aa+1;
            identification_results( n, 6 )=aa; % % the order of non dry period, aa=3 means the third non dry period from the begaining
        elseif identification_results( n, 5 )==0 && identification_results( n-1, 5 )==0
            identification_results( n, 6 )=aa;
        end
    end
end



%% 2. remove
if nargin>4
    identification_results(:,7)=zeros(N,1); % dry or non dry condition after remove minor dry period
    for i=1: max(identification_results( :, 6 ) )
        aa=identification_results( :, 6 )==-i;
        if sum(aa)<REMO 
            identification_results(aa,7)=0;% satisfing remove conditions means it is a minor dry period, don't take into account as drought
        else
            identification_results(aa,7)=1; % 1 means still in dry
        end
    end
    
    % statistics after remove
    identification_results( :,8 )=nan(N,1);
    aa=0;  
    bb=0; 
    
    for n=1:N
        if n==1
            if identification_results( n, 7 )==1
                identification_results( n, 8 )=bb;
            else
                identification_results( n, 8 )=aa;
            end
        else
            if identification_results( n, 7 )==1 && identification_results( n-1, 7 )==1
                identification_results( n, 8 )=bb;
            elseif identification_results( n, 7 )==1 && identification_results( n-1, 7 )==0
                bb=bb-1;
                identification_results( n, 8 )=bb;
            elseif identification_results( n, 7 )==0 && identification_results( n-1, 7 )==1
                aa=aa+1;
                identification_results( n, 8 )=aa;
            elseif identification_results( n, 7 )==0 && identification_results( n-1, 7)==0
                identification_results( n, 8 )=aa;
            end
        end
    end
    
    
    %% 3. Merge
    if nargin>5
        identification_results(:,9)=identification_results(:,7);
        for i=1: max(identification_results( :, 8 ) )
            aa=identification_results( :, 8 )==i; 
            if sum( start_th -  identification_results(aa,4)   ) < MERG
                identification_results(aa,9)=1; 
            else
                identification_results(aa,9)=0; 
            end
        end
        
        % statistics after merge
        identification_results( :, 10 )=nan(N,1);
        aa=0; 
        bb=0; 
        
        for n=1:N
            if n==1
                if identification_results( n, 9 )==1
                    identification_results( n, 10 )=bb;
                else
                    identification_results( n, 10 )=aa;
                end
            else
                if identification_results( n, 9 )==1 && identification_results( n-1, 9 )==1
                    identification_results( n, 10 )=bb;
                elseif identification_results( n, 9 )==1 && identification_results( n-1, 9 )==0
                    bb=bb-1;
                    identification_results( n, 10 )=bb;
                elseif identification_results( n, 9 )==0 && identification_results( n-1, 9 )==1
                    aa=aa+1;
                    identification_results( n, 10 )=aa;
                elseif identification_results( n, 9 )==0 && identification_results( n-1, 9)==0
                    identification_results( n, 10 )=aa;
                end
            end
        end
        
    end
end

%% 4. transfer SPI or SHI value back
identification_results(:,4)=extreme_index0;
% On days when the extreme indicator is NaN, the corresponding identification results should also be set to NaN.
day_is_nan = isnan(extreme_index0);
identification_results(day_is_nan, [5,7,9]) = nan;

end