function [param]=meanMCparam(param)
param.exit=0;    %success!
tolDefault=1e-3;      %absolute error tolerance
nmaxDefault=1e9;       %maximum number of scalar sample points
npcmaxDefault=1e5;       %maximum size of a vector of samples
n0Default=30;      %initial number of samples
n0MinDefault=30;        %minimum allowed initial number of samples
alphaDefault=0.01;      %uncertainty
fudgeDefault=1.1;       %fudge factor for error estimate
timebudgetDefault=30; % default time budget
nbudgetDefault=1e9; % default N budget
%if nargin<3; whfield=''; end
%Error tolerance
if isfield(param,'tol');
    param.tol=abs(param.tol);
else
    param.tol=tolDefault; 
end
% time budget
if isfield(param,'timebudget')
    param.timebudget=abs(param.timebudget);
else
    param.timebudget=timebudgetDefault;
end
%sample budget
if isfield(param,'nbudget')
    param.nbudget=abs(param.nbudget);
else
    param.nbudget=nbudgetDefault;
end

% max N
if isfield(param,'nmax');
    param.nmax=abs(param.nmax);
else
    param.nmax=nmaxDefault;
end
%Maxinum number of scalar values of x per vector
if isfield(param,'npcmax');
    param.npcmax=abs(param.npcmax);
else
    param.npcmax=npcmaxDefault;
end
%Initial sample size
if isfield(param,'n0');
    param.n0=max(param.n0,n0MinDefault);
else
    param.n0=n0Default; 
end
%Fudge factor to multiply error estimate by
if isfield(param,'fudge');
    param.fudge=max(abs(param.fudge));
else
    param.fudge=fudgeDefault; 
end
%Uncertainty
if isfield(param,'alpha');
    param.alpha=min(max(abs(param.alpha),0),1);
else
    param.alpha=alphaDefault; 
end
end