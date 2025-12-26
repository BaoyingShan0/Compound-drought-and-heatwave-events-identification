function [ind, nse, models, criteria_gev, ksTest_best_dist]=best_fitting(data, for_plot, for_ks, for_criteria)


% 23.04.04: 10 to 9 distributions

% aim: find the best distribution for data based on minimum AIC, BIC, RMSD, or/and pass
% ks test.
% test if for_ks=true
%for_criteria: given the criteria for choosing the best distribution, the
%default are AIC, RMSD; 

% input:
% data: a sequence of random variable
% for_plot: if for_plot='p', plot the simulated cdf nd empirical cdf for each distribution; other else, no plot
% for_ks: if for_ks=true, then the best dist should also pass ks test at the
% same time;

% output:
% ind: the index of the best-fitting in models
% nse: the NSE of the best-fitting, -1-1
% models: all the distribution for chosen
% criteria_gev: the MLE, AIC, BIC, NSE for GEV
% ksTest_best_dist: whether the best-fitting pass the KS test

arguments %default value
    data double =[]
    for_plot string = [] %
    for_ks double = [] %
    for_criteria string = ["AIC"; "RMSD"] %
end

%%
% Load the critical value table: the size of table  is 2000*6. The value of each columns is at
% significant 0.001,0.01,0.05,0.1,0.15,0.2, respectively rows mean the sample num from 1 to 2000
% Dn0= readtable('Dn.csv', 'HeaderLines' ,1, 'ReadRowNames' , 1); % skip the first line: the variable name line.
load Dn0.mat
models=[ "Normal","Exponential","Gamma", "gev","InverseGaussian",...
    "logistic", "Loglogistic","Lognormal","Burr", "NOBEST" ];% , "gp"]; % , "P-III"]; %gpÔÚmleÊ±±¨´í
k=[2, 1, 2, 3, 2, 2, 2, 2, 3];% the number of parameter

%%
N=size(models, 2);

MLE=nan(N, 1); NSE=zeros(N,1); AIC=nan(N,1 ); BIC=nan(N,1); h_ks=nan(N,1); RMSD=nan(N,1);
for i=1:N-1
    try % if an error, then skip
        lastwarn('', ''); % reset the lastwarn message and id
        
%         if models(i)=="gp" % for GP
%             data= data( data> min(data) ) - min(data); % gp is the last model, so we can change the value of data which would not influence on other models fitting
%             % could change the data value and then change the MLE
%         end
        
        pd=fitdist(data, models(i) );   % [shape, scale ¦Ò, location ¦Ì]
        % now if a warning was raised, warnMsg and warnId will not be empty.
        [~, warnId] = lastwarn();
        
        % you can check the warning message or id, or just throw the warning as an error if desired
        if(isempty(warnId)) % if there is no warning for fitting, then calculate MLE
            
            MLE(i)=-negloglik(pd); %log likelihood
            AIC(i)=2*k(i)-2*MLE(i);
            BIC(i)=k(i)*log( length(data) ) -2*MLE(i);
            
            
            % for comparing
            % ecdf vs modeled cdf
            [resCDF,queryRes] = ecdf(data );
             % for RMSD, NSE, Nash-Sutcliff efficiency(NSE)
            invcdf=icdf(pd, resCDF);
            
            queryRes=queryRes(2:end-1);  invcdf=invcdf(2:end-1); % no P=1 or P=0;
            %scatter(queryRes, invcdf);
            delta=queryRes-invcdf;
            %  mean square deviation of variable
            RMSD(i) = sqrt( sum(delta.^2)/ length(invcdf) ); %RMSE£º root
            NSE(i)=1- sum(delta.^2)/ sum(  ( invcdf-mean(invcdf) ) .^2) ;

            %% ks test
%             h_ks(i) = ksTestby(data, 0.05, pd ,Dn0)
            
            %% for plotting         
            if for_plot=='p'
                sampleRes = linspace(min(data ),...
                    max(data ),100).';
                fitCDF=cdf(pd, sampleRes );
                figure(i)
                stairs(queryRes,resCDF(2:end-1),'o-')
                hold on
                plot( sampleRes ,fitCDF, 'LineWidth', 1.5)
                xlabel('Variable')
                ylabel('Cumulative probability')
                % title('Empirical and Fitted CDFs');
                title( [ 'log likehood: ' num2str(MLE(i),4) ] )
                grid on
                aa=models(i)+" CDF" ;
                legend({'Empirical CDF', aa},...
                    'Location','southeast')
                hold off;
            end
            
            %

            
        else % skip this distribution
        end
        
    catch
%         i=i+1;
    end
end

% one or more criteria to pick the best distribution
NN =size( for_criteria, 1);
IND=nan(NN,1);
for j =1:NN
    switch for_criteria(j)
        case "AIC"
            [~, IND(j)]=nanmin(AIC);
        case "BIC"
            [~, IND(j) ] =nanmin(BIC);
        case "MLE"
            [~, IND(j)] =nanmax(MLE);
        case "RMSD"
            [~, IND(j) ]=nanmin(RMSD);
    end
end

if  length( unique(IND) )==1 % if the best dist by different criteria is the same
    ind=IND(1);
else
    ind=N; % that means NO BEST
end
    

% % for return
testResults=["not rejected","rejected"];
if ~isnan( h_ks(ind) )
    ksTest_best_dist = testResults( h_ks(ind)+1  ) ;
else
    ksTest_best_dist ="nan";
end
%
nse = NSE(ind);
criteria_gev = [MLE(4), AIC(4), BIC(4), NSE(4) ];
BestModels=models(ind);


%% ks test
function [h, testResults] =ksTestby(data, p, pd ,Dn0)

% 2022.1.11
% this function is to test if data follow give model diatribution by ks test
% at p significant level

% data: one colunm is one series of data
% p: significant level
% model: given distribution model, like "gev", "gpd"

N =size(data,2);

ps=[0.001	0.01	0.05	0.1	0.15	0.2];
Dn=Dn0(:, ps==p);

for i =1:N
    % exclude nan value
    x = data(:,i);
    x(isnan(x))=[];
    
    % get theritical cdf based on given pd
    xx=sort(x) ;
    theop = cdf(pd, xx) ;
    
    % empirical
    n=length(x);
    obspu = [1:1: n]/n;
    obspl =[0:1:n-1] / n;
    for_ks = max( max( abs(theop - obspu'), abs(theop - obspl') ) );
    
    dn=Dn{n,1};
    if for_ks <= dn
        testResults{i}="not rejected"; h(i)=0;
    else
        testResults{i}="rejected"; h(i)=1;
    end
    
end
