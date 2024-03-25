function events = daily_2_events(extreme_daily, extreme_type)

% daily index is gatherd to events:
% events=[num, duration, severity, duration, the duration of the non-dry period i,
% the inter-arrival time, start date of year, month, day, end dates of year, month, day, the year for the middle of the duration]

arguments %default value
    extreme_daily double =[];
    extreme_type string = [];
end



events=nan(max( extreme_daily(:, end) ) ,1) ;
for i=1: -nanmin( extreme_daily(:, end) ) %
    aa=extreme_daily( extreme_daily(:, end )==-i,:) ;
    events(i,1:4)=[i, size(aa,1), sum(aa(:,4) ), sum(aa(:,4) )/size(aa,1)] ; %num, duration, severity, duration
    events(i, 7:12) = [  aa(1, 1:3),  aa(end, 1:3) ];
    if aa(1,1)==aa(end,1)
        events(i, 13) = aa(end,1);
    else
        nn=1;
        for yr= aa(1,1): aa(end,1)
            pp(1,nn) = yr;
            pp(2,nn)=sum( extreme_daily(:, 1 )==yr &  extreme_daily(:, end )==-i  ) ;
            nn=nn+1;
        end
        [~,Idx]=max(pp( 2, : ));
        events(i, 13) = pp(1, Idx);
    end

    if i<2
        events(i,5:6)=nan;
    else
        bb=find(extreme_daily(:, end )==-i  ) ;
        cc=find( extreme_daily(:, end )==-i+1 ); % last events
        events(i,5)=bb(1)-cc(end)-1; % the duration of the non-heatwave period i
        events(i,6)=events(i, 2)+events(i, 5); % the inter-arrival time
    end
end


if extreme_type == "d" | extreme_type =="c" % negative
    events(:,3:4)= - events(:,3:4); % severity and intensity
end


