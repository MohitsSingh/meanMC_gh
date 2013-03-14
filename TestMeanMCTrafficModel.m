<<<<<<< .merge_file_IiNLIV
% this is the drive file to test meanMC package,
% using the Nagel-Schreckenberg traffic model 
param.timebudget=10;% given time budget 
param.nbudget=1e6;% given sample budget
param.n0=70;% given intitial sample szie to estimate sigma
param.fudge=1.1;%given fudge factor
param.alpha=0.01;% given uncertainty
param.tol=1e-2;%given error tolerance
param.npcmax=90;% given piecewise maximum to calculate the mu
=======
% Testing meanMC on the Nagel-Schreckenberg Traffic Model
param.timebudget=20;
param.nbudget=1e6;
param.n0=70;
param.fudge=1.1;
param.alpha=0.01; %uncertainty
param.tol=3e-3; %absolute error tolerance
%param.nmax=1e4;
param.npcmax=90;
>>>>>>> .merge_file_rorNve
[mu,param]=meanMC(@Ytrafficmodel,param)

 