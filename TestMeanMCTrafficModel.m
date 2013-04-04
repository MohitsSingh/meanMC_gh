% this is the drive file to test meanMC package,
% using the Nagel-Schreckenberg traffic model 
param.timebudget=10;% given time budget 
param.nbudget=1e6;% given sample budget
param.n0=70;% given intitial sample szie to estimate sigma
param.fudge=1.1;%given fudge factor
param.alpha=0.01;% given uncertainty
param.tol=1e-2;%given error tolerance
param.npcmax=90;% given piecewise maximum to calculate the mu
[mu,param]=meanMC(@Ytrafficmodel,param)
