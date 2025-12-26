function [h, testResults] =KsTest(data, p, pd , Dn0, Dn0_level)

% 2023.4.03
% this function is to test if data follow give model diatribution by ks test
% at p significant level

% data: one colunm is one series of data
% p: significant level
% model: given distribution model, like "gev", "gpd"

N =size(data,2);
Dn=Dn0(:, Dn0_level==p);

for ii =1:N
    % exclude nan value
    x = data(:,ii);
    x(isnan(x))=[];

    % get theritical cdf based on given pd
    xx=sort(x) ;
    theop = cdf(pd, xx) ;

    % empirical
    n=length(x);
    obspu = [1:1: n]/n;
    obspl =[0:1:n-1] / n;
    ks = max( max( abs(theop - obspu'), abs(theop - obspl') ) );

    dn=Dn(n);
    if ks <= dn
        testResults{ii}="not rejected"; 
        h(ii)=0;
    else
        testResults{ii}="rejected"; 
        h(ii)=1;
    end

end
