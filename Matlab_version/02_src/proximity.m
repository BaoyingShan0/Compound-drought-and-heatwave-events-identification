
function Proximity = proximity(extreme_type, Date, index0, start_th, end_th, REMO, MERG)

% 2023.06.15
% calculate the proximity
% extreme_type == ["d", "p","h","c"] is for drought, pluvial, heatwave and
% coldwave

if extreme_type == "d" || extreme_type =="c"
    index = - index0;
    start_th = -start_th;
elseif extreme_type == "p" || extreme_type =="h"
    index = index0;
end

extreme_daily = PRM_extreme_identification( extreme_type, Date, index0, start_th,end_th ,...
    REMO, MERG);

Proximity=nan(max(extreme_daily( :, end ) ),1);
for i=1: max(extreme_daily( :, end ) )
    aa=extreme_daily( :, end )==i;
    Proximity(i) = sum( start_th -  index(aa,1) );
end
