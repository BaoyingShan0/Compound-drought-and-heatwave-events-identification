function boundarys = bound_remove_merge( Dr_Hw, DATE, index, start_th,end_th, ths, max_annual_days, min_total_number  )

% this function could get a rought boundary thresholds for removing and merging to imporve the computing
% efficiency in main4_robustness

% input:
% max_annual_days: annual days in drought should less than a large number, like 120 days, it's a wide boundary and reasonable limit
% min_total_number: total events number in drought should large than a small number, like 10, it's a wide boundary and reasonable limit

N=length(ths);
years=DATE(end,1)-DATE(1,1)+1;

if Dr_Hw=="Dr"
    min_bound=1;
    drought_daily = PRM2_drought_identification( DATE, index, start_th, end_th,  ths(min_bound), 0 );
    annual_days = sum( drought_daily(:, end-1) )/ years;
    while annual_days > max_annual_days % annual days in drought should less than max_annual_days, it's a wide boundary and reasonable limit
        min_bound=min_bound+1;
        drought_daily = PRM2_drought_identification( DATE, index, start_th, end_th,  ths(min_bound), 0 );
        annual_days = sum( drought_daily(:, end-1) )/ years;
    end
    
    max_bound=N; drought_daily = PRM2_drought_identification( DATE, index, start_th, end_th,  ths(max_bound), 0 );
    total_number = abs(  drought_daily(end, end) );
    while total_number <min_total_number % total events number in drought should large than min_total_number, it's a wide boundary and reasonable limit
        max_bound=max_bound-1;
        drought_daily = PRM2_drought_identification( DATE, index, start_th, end_th,  ths(max_bound), 0 );
        total_number = abs(  drought_daily(end, end) );
    end
elseif Dr_Hw=="Hw"
    min_bound=1;
    heatwave_daily = PRM2_heatwave_identification( DATE, index, start_th, end_th,  ths(min_bound), 0 );
    annual_days = sum( heatwave_daily(:, end-1) )/ years;
    while annual_days > max_annual_days % annual days in drought should less than max_annual_days, it's a wide boundary and reasonable limit
        min_bound=min_bound+1;
        heatwave_daily = PRM2_heatwave_identification( DATE, index, start_th, end_th,  ths(min_bound), -inf );
        annual_days = sum( heatwave_daily(:, end-1) )/ years ;
    end
    
    max_bound=N; heatwave_daily = PRM2_heatwave_identification( DATE, index, start_th, end_th,  ths(max_bound), 0 );
    total_number = abs(  heatwave_daily(end, end) );
    while total_number <min_total_number % total events number in drought should large than min_total_number, it's a wide boundary and reasonable limit
        max_bound=max_bound-1;
        heatwave_daily = PRM2_heatwave_identification( DATE, index, start_th, end_th,  ths(max_bound), -inf );
        total_number = abs(  heatwave_daily(end, end) );
    end
    
end

boundarys =[ths( min_bound ), ths( max_bound) ];

% for i=1:N
%      heatwave_daily = PRM2_heatwave_identification( DATE, index, start_th, end_th,  ths(i), -inf );
%       annual_days(i) = sum( heatwave_daily(:, end-1) )/ years;
%       total_number(i) = abs(  heatwave_daily(end, end) );
% end
% figure; plot(annual_days, total_number)
