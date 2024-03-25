function [best_combination, criteria]=...
    remo_merg(extreme_type, Date, SPI_or_SHI, scale, p, start_th, end_th, events_number_control, events_days_control)

% 2023.06.15
% improve the efficiency, not go throught one by one, only go through the number of
% durations and deficiency area

% 2023.04.07
% works for four extremes

% 2023.04.04
% improve the searching stratage, readue the run time
% this function could get the final thresholds for removing and
% merging

% this is an integer programmming problem with objective of geeting maximum number of extreme events
% under the constraints of below
% 1. severity of events should follow the GEV ditribution
% 2. inter arrival times should follow the exponential distribution
% and some basic requirments
% 3. events number per year should >= events_number_control, like
% 0.1 (10-years return period), otherwise too extreme and no enough data.
% 4. days in events per year should <= events_days_control, like 120 days/yr (3 months), other wise too frequent
% 5. merg_d <= 3*remo_d, Observation summed up

% algarithm for iteration

% input:
% Date: [year, month, day]; 365 days each year
% SPI: daily value, correpondeing with date
% scale:
% p: significant level for ks test
% start_th_d, end_th_d: start and thresholds for dr
% events_number_control, events_days_control

% output:
% [best_combination, criteria]


arguments %default value
    extreme_type string =[]
    Date double =[]
    SPI_or_SHI double =[]
    scale double = []
    p double =0.05
    start_th double = -1
    end_th double = -1
    events_number_control double= 0.1 % events number should >= 0.1 events per year, otherwise too extreme and no eought data
    events_days_control double= 120 % events days should <= 120 days per year, other wise too frequent
end

best_combination=[nan,nan,nan,nan];
criteria="nan";

%%
sig_p=0.05;
load Dn0.mat

%% possible combinations
% works for larger scales
extreme_daily = PRM_extreme_identification( extreme_type, Date, SPI_or_SHI, start_th,end_th ,...
    -inf, -inf );
events = daily_2_events(extreme_daily, extreme_type );
duration = events(:,2);
[ecdf_remo,REMO_d] = ecdf( sort( duration +1 ) );
REMO_d(1)=1; 
ecdf_remo(1)=0;
cc = REMO_d > 5*scale;
REMO_d(cc)=[]; ecdf_remo(cc)=[];
REMO_d=[-inf; REMO_d];
ID_R=length(REMO_d);


%%
REMO_MERG_d=[];  ecdf_remo_merg=[]; LEN=0;
MERG_D= cell(ID_R, 1);
for id_r=1:ID_R
    % id_r
    Proximity = proximity( extreme_type, Date, SPI_or_SHI, start_th,end_th ,...
        REMO_d(id_r), -inf );
    if ~isempty(Proximity)
    %     MERG_d{id_r} = unique( sort( ceil(Proximity) ) )
    [ecdf_merg, merg_d] = ecdf(sort( ceil(Proximity) ));
    merg_d(1) = -inf;
    ecdf_merg(1) = 0;

    cc = merg_d > 3*REMO_d(id_r) | merg_d > 5*scale;
    merg_d(cc) = []; ecdf_merg(cc) = [];
    MERG_D{id_r, 1} = merg_d;
    LEN=max(LEN, length(merg_d));
    end
end

MERG_d=nan(ID_R, LEN);
for id_r=1:ID_R
    cc = MERG_D{id_r, 1};
    MERG_d(id_r, 1:length(cc)) = cc;
end

ID_M=size(MERG_d, 2);
mat_id=nan(ID_R, ID_M); % 0 means impossible, 0.5 mean possible, 1 means best
% give a condition that the merging thresholds are usually small, 
mat_id(isnan( MERG_d ))=0;

%% 
% string matrix for save the the best distribution for fitting severity
Best_dist_dr = strings(ID_R,ID_M);
Best_dist_dr(:) = "NaN";
H_dr = NaN(ID_R,ID_M); % save whether inter arrival time reject the exponential distribution, 0 means no rejection

years=size(Date,1)/365;
criterias={["AIC"; "RMSD"], ["AIC"]}; 
% criterias={["AIC"]};
for ccrieria = 1: size(criterias,2)
    criteria=criterias{ccrieria};

    if ccrieria==2
        mat_id(Best_dist_dr=="NOBEST")=nan;
    end
    % initial value
    % start with [remove_d = scale, rmerge_d = 0]
    id_r= 1; id_m= 1; %
    times=0;
%     tic
    while any(isnan(mat_id(:)))
        % all combinations are searched, 如果还包含nan, keep searching
        times = times + 1; % 1/10 of total 60*60
        [ REMO_d(id_r),  MERG_d(id_r, id_m) ];
        extreme_daily = PRM_extreme_identification( extreme_type, Date, SPI_or_SHI, start_th,end_th ,...
            REMO_d(id_r), MERG_d(id_r, id_m) );
%         heatwave_daily0=PRM_extreme_identification( "hw", Date, SHI, start_th_h, end_th_h, best_combination_hw(1), best_combination_hw(2));

        % mat_id(id_r, id_m )
        events_per_yr = -min( extreme_daily(:,end)) / years;
        days_per_yr= sum( extreme_daily(:,end-1)  ) / years;


        % events_number_control, events_days_control
        if events_per_yr < events_number_control % if the events number are already smaller than the events_number_control,
            % then the combinations with the same merging threhold and large remove
            % threshold will have fewer number of events, 一列下方的所有combinations are
            % impossible
            mat_id(id_r:end , id_m) =0;
        elseif days_per_yr > events_days_control % if the days per year are already larger than the events_days_control
            mat_id(id_r, id_m:end ) = 0; % then the combinations with the same removing threhold and larger merging threshold all have the more events days,
            % then [id_r, id_m ] 右侧一行的所有combinations are not possible
        else
            [best_dist_dr, h_dr] = meet_two_assumption_or_not( extreme_type, extreme_daily, criteria, p, Dn0, Dn0_level);
            Best_dist_dr(id_r, id_m) =best_dist_dr;
            H_dr(id_r, id_m)=h_dr;

            if best_dist_dr=="gev" && h_dr==0  % if best_dist_dr=="gev" & h_dr==0, then meet the two assumptions
                mat_id(id_r:end, id_m:end ) = 0;  % as the maximum events number is the final aim
                % for the events number [7 1] > [8 1] > [9 1], and [7 1] > [7 2] >
                % [7 3], so the [id_r, id_m] All combinations in the matrix box in the lower right corner cannot have [id_r, id_m], and all combinations in the rectangular box are 0
                mat_id(id_r, id_m ) = 0.5;  % this combination is a possible one
            else
                mat_id(id_r, id_m ) = 0;
            end
        end


        % main question, how to decide the direction for testing the next combination?
        % Fixed the same merg_d, first search in the direction of small remov_d, if not, then search in the direction of large remov_d;
        % After the merge_d ends, go to merg_d+1
        % Definitely not efficient enough
        while id_r<= ID_R && id_m<= ID_M && ~isnan( mat_id(id_r, id_m ) )
            if any(isnan( mat_id(:, id_m ) )) % if the same merge_d is still searched
                for_search= isnan( mat_id(:, id_m ) );
                if any( for_search(1:id_r) ) 
                    pp=find(for_search(1:id_r)==1);
                    id_r=pp(end);
                elseif any( for_search(id_r:end) )
                    for_search(1:id_r)=0;
                    pp=find(for_search==1);
                    id_r=pp(1);
                end
            else % else The same merg_d has been searched, go to merg_d+1(id_m +1)
                id_m=id_m+1;
            end
        end
    end
%     toc

    % pick the best according to the max number of events
    [id_r_best,id_m_best] = find( Best_dist_dr=="gev" & H_dr==0);

    if ~isempty(id_r_best) % if there is no best using the AIC and RMSD as the criteria, then only go to the AIC
        N0=length(id_r_best);

        for i=1:N0
            extreme_daily = PRM_extreme_identification( extreme_type, Date, SPI_or_SHI, start_th,end_th ,...
            REMO_d(id_r_best(i)), MERG_d(id_r_best, id_m_best(i)) );
            events_per_yr(i) = -min( extreme_daily(:,end)) / years;
            days_per_yr(i) = sum( extreme_daily(:,end-1)  ) / years;
        end

        [~, ind]=max(events_per_yr);
        best_combination=[REMO_d(id_r_best(ind))', MERG_d(id_r_best(ind), id_m_best(ind))' ,events_per_yr(ind)',days_per_yr(ind)'];

        break
    end
end

