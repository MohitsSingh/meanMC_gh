function mu=gaussianfunDoctest(a)
% test gaussian fun
%
% gaussianfunDoctest(a)
%    returns (mu)
%
% Examples:
%
% >> gaussianfunDoctest(1)
% 
% ans =
% 
%     1
% 
% >> gaussianfunDoctest(100)
% 
% ans =
% 
%     1.0**
% 
% >> gaussianfunDoctest('hi')
% ??? Error using ==> gaussianfunDoctest ***
% gaussianfunDoctest(a) requires value to be a number
% 
%
% TWO blank lines before the prose description of the function continues
%

format long
param.timebudget=100;% given time budget 
param.nbudget=1e6;% given sample budget
param.n0=70;% given intitial sample szie to estimate sigma
param.fudge=1.1;%given fudge factor
param.alpha=0.01;% given uncertainty
param.tol=1e-2;%given error tolerance
param.npcmax=90;% given piecewise maximum to calculate the mu
y=@(n) exp(-a^2.*(rand(n,1)-pi/10).^2)./(sqrt(pi)/(2*a)*(erf(a-pi*a/10)-erf(-pi*a/10)));
mu=meanMC(@(n)y(n),param);
mu=round(50*mu)/50;
end


 
