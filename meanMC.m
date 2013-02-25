function [mu,param]=meanMC(Yinput,param)
tstart=tic;
[param]=meanMCparam(param); %check validity of inputs
%% convert time budget to sample budget
r_data=[2,10]; % intital replication number to calculate
polydegree=1; % set the polyfit degree 
[param.ntpredict,~,param.timebudget,~]=...
    timeBudgetToSampleBudget(param.timebudget,r_data,polydegree);
% use the above function to convert the time budget to sample budget
function [r_predict,t_data,t_budget,p]=timeBudgetToSampleBudget...
        (t_budget,r_data,polydegree)
r1=1;
tic;
Ytrafficmodel(r1);
t1=toc;  % let it run once to load all the data.warm up the machine.
clear t1;

N = length(r_data); % number of repilcations used to estimate
t_data=zeros(size(r_data)); 

for j=1:N 
tic
Ytrafficmodel(r_data(j));
t_data(j)=toc; % calculate the time we used when given different replication number
end
p=polyfit(r_data,t_data,polydegree); % using polyfit to fit all the r and t data
r_predict=ceil(fzero(@(R) polyval(p,R)-t_budget,1e-6)); % predict the replication
% number we need to use when given a time budget
end
%% main part of MeanMC 
param.nmax=min(param.ntpredict,param.nbudget);% set the minimum n budget
param.nmumax=max((param.nmax-param.n0),1);% get the n budget to estimate mean
if param.nmumax < param.n0;%if the cost to estimate the variance is too big, just
    param.exit=2; %estimate the mean,skip stage one.
    meanMCerr(param,tstart);%print warning message
    mu=mean(Yinput(param.nmumax+param.n0)); 
else
    Yval=Yinput(param.n0); %use the initial sample size n0 to estimate variance
    param.var=var(Yval);% calculate the sample variance--stage 1
    sig0=sqrt(param.var);
    sig0up=param.fudge*sig0;
    alpha1=1-sqrt(1-param.alpha);
    param.kurtmax=(param.n0-3)/(param.n0-1) ...
        + ((alpha1*param.n0)/(1-alpha1))*(1-1/param.fudge^2)^2;%get kmax
    if sig0up==0; % if the variance is zero, just take one sample 
        param.n=param.n0;
    else
        toloversig=param.tol/sig0up;%tolerance over sigma
        ncheb=ceil(1/(alpha1*toloversig.^2)); % use chebyshev inequality to estimate n
        A=18.1139;
        A1=0.3328;
        A2=0.429; %three constants in Berry-Esseen inequality
        M3upper=param.kurtmax^(3/4);%using jensen inequality to 
                                %bound the third moment
        BEfun=@(logsqrtn)normcdf(-exp(logsqrtn).*toloversig)...
            +exp(-logsqrtn).*min(A1*(M3upper+A2), ...% Berry-Esseen Inequality
            A*M3upper./(1+(exp(logsqrtn).*toloversig).^3))- alpha1/2; 
        if BEfun(log(sqrt(param.n0))) <= 0 || ncheb <=param.n0;%n0 is too large
            param.n=param.n0;
        else
        logsqrtnCLT=log(norminv(1-alpha1/2)/toloversig);
        param.n=min(ncheb,ceil(exp(2*fzero(BEfun,logsqrtnCLT))));
        % get the min n by using cheb and BEfun
        end
    end
    if param.n > param.nmumax  %too many samples for the time allowed
    param.exit=1;
    meanMCerr(param,tstart); %print warning message
    param.n=floor(param.nmumax);
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
    param.mu=sumY/param.n; %approximate integral is sample mean
    param.n=param.n+param.n0; %total number of samples used
    mu=param.mu; %assign answer
end
param.time=toc(tstart); %elapsed time
end