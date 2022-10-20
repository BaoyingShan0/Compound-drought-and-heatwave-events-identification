function compound =identify_comound (Index1_0, Index2_0, Index1, Index2,  type)

%2022.7.26
% identify compound events of two extremes: extreme1, extreme2

% inputs: 
% Index1_0, index2_0:  are nan for non extreme days of two extremes
%  Index1, Index2: daily index values of two extremes
% types=1,2,3,4 is for union, conditions on extreme 1, condition on extreme 2, and intersecion.

% output:
% compound: the first column is the logical value of this in compound
% events or not; the second column is the order of compound events


N=size(Index1_0, 1);
idx_1=nan(N, 1);  types=["Union", "ex1|ex2", "ex2|ex1", "intersection"];

switch type
    case 1 % union
        for i=1:N
            if isnan( Index1_0(i) ) & isnan(Index2_0(i) )
                idx_1(i, 1)=0;
            else
                idx_1(i, 1)=1;
            end, end
        
    case 2 %extreme1 as condition
        for i=1:N
            if isnan(Index1_0(i) )
                idx_1(i, 1)=0;
            else
                idx_1(i, 1)=1;
            end, end
        
    case 3 % extreme2 as condition
        for i=1:N
            if isnan(Index2_0(i) )
                idx_1(i, 1)=0;
            else
                idx_1(i, 1)=1;
            end, end
        
    case 4 % ingtersection
        for i=1:N
            if isnan(Index1_0(i) ) ||  isnan(Index2_0(i) )
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

%% to idx, % ignoring the univariate events 
idx(:, 1:2)=zeros(N,2) ;
% ignoring the univariate events
for i=1: max( idx_1(:, 2) ) 
    aa= find( idx_1(:, 2)==i ) ;
    if  ~all( isnan(  Index1_0(aa) -Index2_0(aa)  ) ) % ignoring the univariate events % nansum(SPI_0(aa) ) ~=0 &  nansum(SHI_0(aa) ) ~=0 ;
        idx(aa,1)=1;
        idx(aa, 2)=max( idx(:, 2) )+1; % update number of compound events
    end 
end

% %% final compound events
% a1=size(idx, 2);
% compound_events=[] ; idx(:, a1+1:a1+2)=0;
% for i=1: max( idx(:, a1) ) %
%     aa= find( idx(:, a1)==i ) ;
%     
%     if  ~all( isnan(  Index1_0(aa) -Index2_0(aa)  ) ) % ignoring the univariate events % nansum(SPI_0(aa) ) ~=0 &  nansum(SHI_0(aa) ) ~=0 ;
%         if size(compound_events,1)==0
%             intervals=nan;
%         else
%             tt=size(compound_events,1);
%             intervals=aa(1)-compound_events(tt, 6)+1 ;
%         end
%         
%         %key part: sum all values in duration for marginal severity
%         severity1=  nansum( Index1(aa) ); severity2=  nansum( Index2(aa) );
%         if nanmean(Index1_0)<0
%             severity1=  -nansum( Index1(aa) );
%         end
%         if nanmean(Index2_0)<0
%             severity2=  -nansum( Index2(aa) );
%         end
%         
%         compound_events=[compound_events; length(aa), severity1,  severity2, intervals, aa(1), aa(end)] ; % duration, severity of drought, severity of heatwave, intervals, 开始和结束的序号，这个序号对应着日期
% 
%         idx(aa, a1+1)=1;
%         idx(aa, a1+2)=max( idx(:, a1+2) )+1; % update number of compound events
%     end
%     
% end

compound=idx(:, end-1:end); % for return

