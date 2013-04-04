function [mu,param]=meanMC(Yinput,param)
%  meanMC evaluate the mean of the given data to a specific error tolerence.

%  param.n0: Initial given sample size to estimate the sample variance
%  param.nbudget: the maximum sample budget to do the two-stage estimation
%  param.timebudget: the maximum time budget to do the two-stage estimation
%  param.ntpredict:Using the given time budget(param.timebudget) to ...
%                  estimate how many samples(param.ntpredict) we could afford
%                  to do the two-stage estimation
%  param.nmax:take the minimum of param.nbudget and param.ntpredict
%              as our sample budget to do the two-stage esimation
%             
%  param.ntpredict2:adding another data from n0 replication and do the
%                    estimation again to get how many sample budget left 
%                    to estimate the mu(second-stage estimation)
%  param.nmumax:the maximum sample budget to estimate mu (second-stage estimation)
%  Yinput:the given data whos mean we want to estimate

tstart=tic; % start the clock
param=meanMCparam(param); %check validity of inputs
%% convert time budget to sample budget
r_data=[2,10]; % do a small number of replication to train the model
[param.ntpredict,t_data]=timeBudgetToSampleBudget(param.timebudget,r_data,Yinput);
% using a function below to convert the time budget to sample budget

%% main part of MeanMC 
param.nmax=max(max(param.ntpredict,1),param.nbudget);% set the minimum n budget
if  param.n0>param.nmax/2;%if the cost to estimate the variance is too big, just
    param.exit=2; %estimate the mean,skip stage one.
    meanMCerr(param,tstart);%print warning message
    mu=mean(Yinput(param.nmax)); %  calculate the mean
else
    tic
    Yval=Yinput(param.n0); %use the initial sample size n0 to estimate variance
    n0t_data=toc; % get the time used to evaluate n0 samples
    param.var=var(Yval);% calculate the sample variance--stage 1
    sig0=sqrt(param.var);
    sig0up=param.fudge*sig0;
    alpha1=1-sqrt(1-param.alpha);
    param.kurtmax=(param.n0-3)/(param.n0-1) ...
        + ((alpha1*param.n0)/(1-alpha1))*(1-1/param.fudge^2)^2;%get kmax
    if sig0up==0; % if the variance is zero, just take n0 samples 
        param.n=param.n0;
    else
        toloversig=param.tol/sig0up;%tolerance over sigma
        ncheb=ceil(1/(alpha1*toloversig.^2)); % use Chebyshev inequality to estimate n
        A=18.1139;
        A1=0.3328;
        A2=0.429; %three constants in Berry-Esseen inequality
        M3upper=param.kurtmax^(3/4);%using Jensen inequality to 
                                %bound the third moment
        BEfun=@(logsqrtn)normcdf(-exp(logsqrtn).*toloversig)...
            +exp(-logsqrtn).*min(A1*(M3upper+A2), ...% Berry-Esseen Inequality
            A*M3upper./(1+(exp(logsqrtn).*toloversig).^3))- alpha1/2; 
        if BEfun(log(sqrt(param.n0))) <= 0 || ncheb <=param.n0;
            param.n=param.n0;%the Chebyshev n or the BE n is too small, just
                             %use param.n0;
        else
        logsqrtnCLT=log(norminv(1-alpha1/2)/toloversig);
        param.n=min(ncheb,ceil(exp(2*fzero(BEfun,logsqrtnCLT))));
        % get the min n (used to estimate mu) by using cheb and BEfun
        end
    end
    r_data2=[r_data,param.n0];%add another param.n0 data and estimate how many samples
                            % we could afford to estimate the mu--stage 2
    t_data2=[t_data,n0t_data];
    [~,m,b]=regression(r_data2,t_data2); % using linear regression
                                        % to fit all the r and t data
    param.ntpredict2=floor((param.timebudget-toc(tstart)-b)/m); 
    % to estimate sample budget left
    param.nmumax=max(min(param.nmax-param.n0,param.ntpredict2),1);
    % take the mimimum sample budget
    if param.n > param.nmumax % too many samples for the time allowed
    param.exit=1;
    meanMCerr(param,tstart); % print warning message
    param.n=param.nmumax;
    end
%%  Split The Param.n into columns 
    nopt=min(param.npcmax,param.n);
    %Put all samples in one or more vectors
    nn=floor(param.n/nopt); %number of loop steps
    nremain=param.n-nn*nopt; %# of samples in last loop step
    nloop=repmat(nopt,1,nn); %vector of # of samples per loop step
    if nremain>0; nloop=[nloop nremain]; nn=nn+1; end
    sumY=0;
    for iloop=1:nn %loops to save memory  
        sumY=sumY+sum(Yinput(nloop(iloop)));
    end
%%  Estimate mu 
    param.mu=sumY/param.n; %calculate the mean
    param.n=param.n+param.n0; %total number of samples used 
    mu=param.mu; %assign answer
end
param.time=toc(tstart); %elapsed time
end
% function to convert time budget to sample budget
function [r_predict,t_data]=timeBudgetToSampleBudget(t_budget,r_data,Yinput)
r1=1;
tic;
Yinput(r1);
t1=toc;  % let it run once to load all the data.warm up the machine.
clear t1;

N = length(r_data); % number of repilcations used to estimate
t_data=zeros(size(r_data)); 

for j=1:N 
tic
Yinput(r_data(j));
t_data(j)=toc; % calculate the time we used when given different replication numbers
end
[~,M,B]=regression(r_data,t_data); % using linear regression to fit all the r and t data
r_predict=floor((t_budget-B)/M); % predict the sample budget
end