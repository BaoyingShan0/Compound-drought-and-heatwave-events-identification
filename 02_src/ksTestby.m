function [h, testResults] =ksTestby(data, p, model)

% 2022.1.11
% 2022.1.19 fix the error: skip the first indexing column when load Dn0
% this function is to test if data follow give model diatribution by ks test
% at p significant level

% data: one colunm is one series of data
% p: significant level
% model: given distribution model, like "gev", "gpd"

N =size(data,2);

% Load the critical value table: this table 
% size is 2000*6 . the value of each columns is at
% significant 0.001,0.01,0.05,0.1,0.15,0.2, respectively
% rows mean the sample num from 1 to 2000
Dn0= readtable('Dn.csv', 'HeaderLines' ,1, 'ReadRowNames' , 1); % skip the first line: the variable name line. and first column, important
ps=[0.001	0.01	0.05	0.1	0.15	0.2];
Dn=Dn0(:, ps==p);

for i =1:N
    % exclude nan value
    x = data(:,i);
    x(isnan(x))=[];
    
    
    % fit distribution and get theritical cdf
    pd=fitdist(x, model ) ;
    xx=sort(x) ;
    theop = cdf(pd, xx) ;
    
    % empirical
    n=length(x);
    obspu = [1:1: n]/n;
    obspl =[0:1:n-1] / n;
    ks = max( max( abs(theop - obspu'), abs(theop - obspl') ) );
    
    dn=Dn{n,1};
    if ks <= dn
        testResults{i}="not rejected"; h(i)=0;
    else
        testResults{i}="rejected"; h(i)=1;
    end
    
end


