function mu=trafficmodelDoctest(tol)
% test traffic model
%
% trafficmodelDoctest(tol)
%    returns (mu)
%
% Examples:
%
% >> trafficmodelDoctest(1e2)
% 
% ans =
% 
%     4.4***
% 
% >> trafficmodelDoctest(1e3)
% 
% ans =
% 
%     4.45***
% 
% >> trafficmodelDoctest('hi')
% ??? Error using ==> trafficmodelDoctest ***
% trafficmodelDoctest(tol) requires value to be a number
% 
%
% TWO blank lines before the prose description of the function continues
%


param.timebudget=100;% given time budget 
param.nbudget=1e6;% given sample budget
param.n0=70;% given intitial sample szie to estimate sigma
param.fudge=1.1;%given fudge factor
param.alpha=0.01;% given uncertainty
param.tol=tol;%given error tolerance
param.npcmax=90;% given piecewise maximum to calculate the mu
[mu,~]=meanMC(@Ytrafficmodel,param);
end
