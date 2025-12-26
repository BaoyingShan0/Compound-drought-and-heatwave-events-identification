function [compound,compound_daily] =identify_compound (ex1_daily, ex2_daily,  type)

%2022.7.26
% identify compound events of two extremes: extreme1, extreme2
% 2022.12.15
% update the order of four types

% 2023.4.19
% change the input, delete: SPI0 SHI0
% add out put: compound_daily

% inputs: 
% Index1_0, index2_0:  are nan or 0 for non extreme days of two extremes,
% are 1 for extreme days
% Index1, Index2: daily index values of two extremes
% types=["ex1-and-ex2", "ex1-or-ex2", "ex1-cond-ex2", "ex2-cond-ex1"];

% outputs:
% compound: the first column is the logical value of this in compound
% events or not; the second column is the order of compound events
% compound_daily: [Date, ex1_daily during the compound events, ex2_daily during the compound events, compound]


N=size(ex1_daily, 1);
idx_1=nan(N, 1);  types=["ex1-and-ex2", "ex1-or-ex2", "ex1-cond-ex2", "ex2-cond-ex1"];

Index1_0= ex1_daily(:,4); Index1_0(ex1_daily(:, end-1)==0 )=nan; 
Index2_0= ex2_daily(:,4); Index2_0(ex2_daily(:, end-1)==0 )=nan; 

%% 
switch type
        case 1 % ingtersection, ex1-and-ex2
        for i=1:N
            if ( isnan(Index1_0(i) ) ||  isnan(Index2_0(i) ) ) || ( (Index1_0(i)==0 ) ||  (Index2_0(i)==0 ))
                idx_1(i, 1)=0;
            else
                idx_1(i, 1)=1;
            end, end
        
    case 2 % union, ex1-or-ex2
        for i=1:N
            if ( isnan( Index1_0(i) ) & isnan(Index2_0(i) ) ) ||  ( Index1_0(i) ==0  & Index2_0(i) ==0 )
                idx_1(i, 1)=0;
            else
                idx_1(i, 1)=1;
            end, end
        
    case 3 % ex1-cond-ex2
        for i=1:N
            if ( isnan(Index1_0(i) ) ) || ( Index1_0(i) ==0)
                idx_1(i, 1)=0;
            else
                idx_1(i, 1)=1;
            end, end
        
    case 4 % ex2-cond-ex1
        for i=1:N
            if ( isnan(Index2_0(i) ) ) ||  ( Index2_0(i) ==0 )
                idx_1(i, 1)=0;
            else
                idx_1(i, 1)=1;
            end, end
        

end

%% to compound events
aa=0; % the order of non compound period
bb=0; % the order of compound period

for n=1:N
    if n==1
        if idx_1( n, 1 )==1
            idx_1( n, 2 )=bb;
        elseif idx_1( n, 1 )==0
            idx_1( n, 2 )=aa;
        end
    else
        
        if idx_1( n, 1 )==1 && idx_1( n-1, 1 )==1
            idx_1( n, 2 )=bb;
        elseif idx_1( n, 1 )==1 && idx_1( n-1, 1 )==0
            bb=bb+1;
            idx_1( n, 2 )=bb;  % the order of dry period, bb=-2 means the second dry period from the begaining
        elseif idx_1( n, 1 )==0 && idx_1( n-1, 1 )==1
            aa=aa-1;
            idx_1( n, 2 )=aa; % % the order of non dry period, aa=3 means the third non dry period from the begaining
        elseif idx_1( n, 1 )==0 && idx_1( n-1, 1 )==0
            idx_1( n, 2 )=aa;
        end
    end
end

%% compound evetns;
% to idx, % ignoring the univariate events 
idx(:, 1:2)=zeros(N,2) ;
% ignoring the univariate events
for i=1: max( idx_1(:, 2) ) 
    aa= find( idx_1(:, 2)==i ) ;
    if all(~isnan(Index1_0)) & abs( mean (  Index1_0(aa) - Index2_0(aa)  ) ) < 1  
    % ignoring the univariate events % nansum(SPI_0(aa) ) ~=0 &  nansum(SHI_0(aa) ) ~=0 ;
        idx(aa,1)=1;
        idx(aa, 2)=max( idx(:, 2) )+1; % update number of compound events

    elseif  ~all(~isnan(Index1_0)) & ~all( isnan( (  Index1_0(aa) - Index2_0(aa)  ) ) )
        idx(aa,1)=1;
        idx(aa, 2)=max( idx(:, 2) )+1; % update number of compound events
    end
end

% for return
compound=idx(:, end-1:end); % for return, the real compound events

aa1=ex1_daily; aa2 = ex2_daily;
aa1(compound(:,1)==0,:)=0; aa2(compound(:,1)==0,:)=0;
compound_daily=[ex1_daily(:,1:3), aa1(:,4), aa2(:,4), compound];

