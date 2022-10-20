function [idx]=PRM2_drought_identification(Date, index, start_th, end_th, REMO, MERG)


%%
% change 11.3, try the one severity threshold for merging

% Use pre identification, removing, merging three steps to identify drought

%output: idx: each column is [ year, month, day, drought index, dry or non dry
% symbol after pre identification (1 means in dry period and 0 is in non
% dry) , period orders after preidentification (negative is the order of
% dry period and positive is for non-dry period), dry or non dry symbol
% after removing, period orders after removing, dry or non dry symbol
% after merging, period orders after merging]

% input:
%Date: time: 3 columns£¬[year, month, day];
% index: drought index value, tie series;
% STEP 1: pre identification to get dry period and non dry period. a dry period starts when index is less or equal than the start_th
% and end when the index becomes larger than end_th;
%start_th and end_th can be the same value.

% STEP 2: remove or delete minor dry period
% condition 1: If the duration of the dry period is samller than REMO, then
% remove it

% STEP 3: merge the non-independent near dry periods
% if the sum of (SPI+start_th (1) ) during intevals is less then  MERG,
% then merge the two near dry-spells two one
%


%% 1. pre identify

N=size(index,1 );

idx=[Date, ones(N, 3-size(Date,2 ) ), index, nan(N,1) ];

for n=1:N
    if n==1
        if idx(n,4)<=start_th,  idx( n, 5 )=1;
        else, idx( n, 5 )=0; end
    else
        if idx(n,4)<=start_th
            idx( n, 5 )=1; % 1 means in dry period
        elseif idx( n-1, 5 )==1 &&  idx(n,4)<end_th
            idx( n, 5 )=1; % 1 means in dry period
        else
            idx( n, 5 )=0; %0 means in non dry period
        end, end, end

idx( :,6 )=nan(N,1);
aa=0; % the order of non dry period
bb=0; % the order of dry period

% the number of dry periods and wet periods
% all the days in the same period have same number
for n=1:N
    if n==1
        if idx( n, 5 )==1
            idx( n, 6 )=bb;
        elseif idx( n, 5 )==0
            idx( n, 6 )=aa;
        end
    else
        if idx( n, 5 )==1 && idx( n-1, 5 )==1
            idx( n, 6 )=bb;
        elseif idx( n, 5 )==1 && idx( n-1, 5 )==0
            bb=bb-1;
            idx( n, 6 )=bb;  % the order of dry period, bb=-2 means the second dry period from the begaining
        elseif idx( n, 5 )==0 && idx( n-1, 5 )==1
            aa=aa+1;
            idx( n, 6 )=aa; % % the order of non dry period, aa=3 means the third non dry period from the begaining
        elseif idx( n, 5 )==0 && idx( n-1, 5 )==0
            idx( n, 6 )=aa;
        end
    end
end

% cc=nan(N, 1);cc( idx(:,5)==1)=idx( idx(:,5)==1, 4 ); %Fig_Size; bar(t1, cc )
pre_iden_number=[min( idx(:, 6 )) ,max(idx( :, 6 ) ) ]; % total number of, 170 dry period

% cdf of the number of days for dry cells and wet cells
for i=1: max(idx( :, 6 ) )
    wet_cell(i,1)=sum(idx( :, 6 )==i) ;
end
% figure(6)
% ecdf(dry_cell); title('non-dry period'); xlim([0, 100]);
wet_prc_50=prctile(wet_cell,50);

for i=1: max(idx( :, 6 ) )
    dry_cell(i,1)=sum(idx( :, 6 )==-i) ;
end
dry_prc_50=prctile(dry_cell,50);


%% 2. remove
if nargin>4
    idx(:,7)=[zeros(N,1)]; % dry or non dry condition after remove minor dry period
    for i=1: max(idx( :, 6 ) )
        aa=idx( :, 6 )==-i; % find all dry periods
        if sum(aa) < REMO 
            idx(aa,7)=0;% satisfing remove conditions means it is a minor dry period, don't take into account as drought
        else
            idx(aa,7)=1; % 1 means still in dry
        end
    end
    
    % statistics after remove
    idx( :,8 )=nan(N,1);
    aa=0;  % wet period orders
    bb=0; % dry period orders
    
    for n=1:N
        if n==1
            if idx( n, 7 )==1
                idx( n, 8 )=bb;
            else
                idx( n, 8 )=aa;
            end
        else
            if idx( n, 7 )==1 && idx( n-1, 7 )==1
                idx( n, 8 )=bb;
            elseif idx( n, 7 )==1 && idx( n-1, 7 )==0
                bb=bb-1;
                idx( n, 8 )=bb;
            elseif idx( n, 7 )==0 && idx( n-1, 7 )==1
                aa=aa+1;
                idx( n, 8 )=aa;
            elseif idx( n, 7 )==0 && idx( n-1, 7)==0
                idx( n, 8 )=aa;
            end
        end
    end
    
    cc=nan(N, 1); cc( idx(:,7)==1)=idx( idx(:,7)==1, 4 ); % bar(t1, cc )
    after_remove_num=[min( idx(:, 8 )) , max(idx( :, 8 ) ) ] ;% total number of dry period after removing minor dry periods
    
    
    %% 3. Merge
    if nargin>5
        idx(:,9)=idx(:,7);
        for i=1: max(idx( :, 8 ) )
            aa=idx( :, 8 )==i; % find all wet period
            if sum( idx(aa,4) - start_th)<MERG
                idx(aa,9)=1; % this means this small non dry period belongs to the big dry period with before and after dry periods
            else
                idx(aa,9)=0; % still in wet
            end
        end
        
        % statistics after merge
        idx( :, 10 )=nan(N,1);
        aa=0;  %positive: wet period numbers
        bb=0; % negative: dry period numbers
        
        for n=1:N
            if n==1
                if idx( n, 9 )==1
                    idx( n, 10 )=bb;
                else
                    idx( n, 10 )=aa;
                end
            else
                if idx( n, 9 )==1 && idx( n-1, 9 )==1
                    idx( n, 10 )=bb;
                elseif idx( n, 9 )==1 && idx( n-1, 9 )==0
                    bb=bb-1;
                    idx( n, 10 )=bb;
                elseif idx( n, 9 )==0 && idx( n-1, 9 )==1
                    aa=aa+1;
                    idx( n, 10 )=aa;
                elseif idx( n, 9 )==0 && idx( n-1, 9)==0
                    idx( n, 10 )=aa;
                end
            end
        end
        
        after_merge_num=[min( idx(:, 10 )) ,max(idx( :, 10 ) ) ] ;% total number of dry period after merging
        
    end
end

%% Plot the final identified drought events
% dd=size(idx, 2);
% t1=datetime(Date);  cc=nan(N, 1); cc( idx(:,dd-1)==1)=idx( idx(:,dd-1)==1, 4 );
% 
% Fig_Size;plot(t1, idx(:, 4));  datetick('x','yyyy');
% hold on;  plot( t1, start_th*ones(length(idx(:,4) ),1) ,'b--' );
% plot( t1, end_th*ones(length(idx(:,4) ),1)  ,'k--' );
% bar(t1, cc ); %xlim([t1(366), t1(20000)]);
% 
% title(' Final identified drought events');
% xlabel('Date'); ylabel('Drought index');

end