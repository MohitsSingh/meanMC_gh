% Testing meanMC on the Nagel-Schreckenberg Traffic Model
param.timebudget=20;
param.nbudget=1e6;
param.n0=70;
param.fudge=1.1;
param.alpha=0.01; %uncertainty
param.tol=3e-3; %absolute error tolerance
%param.nmax=1e4;
param.npcmax=90;
[mu,param]=meanMC(@Ytrafficmodel,param)

 