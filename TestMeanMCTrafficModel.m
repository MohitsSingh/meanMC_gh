
param.timebudget=1;
param.nbudget=1e6;
param.n0=70;
param.fudge=1.1;
param.alpha=0.01;
param.tol=1e-2;
%param.nmax=1e4;
param.npcmax=90;
[mu,param]=meanMC(@Ytrafficmodel,param)

 