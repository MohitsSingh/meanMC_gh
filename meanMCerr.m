function [param,mu]=meanMCerr(param,tstart)
%Handles errors in meanMC and meanMCparam
%to give an exit with information
%param.exit = 0   success
%             1   too many samples required
%             2   initial sample size is too large
if ~isfield(param,'exit'); return; end
if param.exit==0; return; end
switch param.exit
    case 1 %too many samples
        fprintf(2,['Warning: tried to evaluate at ' int2str(param.n) ...
            ' samples,\n' ...
            'which is more than the allowed maximum of ' ...
            num2str(param.nmax) ' samples\n']);
            return
    case 2 % nmax too small, n0 exceeds half of the total sample budget
        fprintf(2,['Warning: n0 = ' int2str(param.n0) ...
            ' which exceeds half of '...
            'nmax = ' int2str(param.nmax) ...
            ',\n just use nmax = ' int2str(param.nmax)...
            ' samples to estimate mean,'...
            ' skip stage 1']);
           return
end
param.mu=NaN;
mu=param.mu;
if nargin>1; param.time=toc(tstart); end

