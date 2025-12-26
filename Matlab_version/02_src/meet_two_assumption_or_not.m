function  [best_dist, h] =meet_two_assumption_or_not(extreme_type, extreme_daily, criteria, p, Dn0, Dn0_level)


M=  -min( extreme_daily(:,end));
events=nan(M,5);
for i=1:M %
    aa=extreme_daily( extreme_daily(:, end )==- i,:) ;
    events(i,1:4)=[i, size(aa,1), sum(aa(:,4) ), sum(aa(:,4) )/size(aa,1)] ; %if extreme_type == "wet" || "hw" % positive
    %num, duration, severity, intensity
    if i<2
        events(i,5)=nan;
    else
        bb=find(extreme_daily(:, end )==-i  ) ;
        cc=find( extreme_daily(:, end )==-i+1 ); % last events
        events(i,5)=bb(1)-cc(end)-1; % the duration of the non-dry period i
    end
end
arrivals_spells=events(2:end, 2)+events(2:end, 5);

if extreme_type == "d" || extreme_type =="c" % negative
    events(:,3:4) = -events(:,3:4);% severity and intesnsity
end

% Assumption 1：arrival should follow Poisson process, i.e. inter
% arrival time should follow exponential distribution, us ks test
pd = fitdist(arrivals_spells, "exp" ) ;
h = KsTest(arrivals_spells, p, pd, Dn0, Dn0_level); % if h==0,  arrivals follow exponential
best_dist="nan";
if h==0
    % Assumption 2：severity should follow GEV
    % is the best distribution for severity GEV? based on AIC and RMSD
    pd = fitdist(events(:, 3), "gev" );
    h2 = KsTest(events(:, 3), p, pd, Dn0, Dn0_level); % if h_2==0, severity follow exponential
    if h2==0
%         best_dist="gev"; % try not comparing with other distributions，
%         perform nice for heatwave, but not for drought, with accumulation
%         period 45 days, 60 days, 90 days, very small removing
%         threhsolds, like 1, 2, 3 
        [ind, ~, models]=best_fitting( events(:, 3), "", true, criteria );
        best_dist=  models( ind );
    end
end