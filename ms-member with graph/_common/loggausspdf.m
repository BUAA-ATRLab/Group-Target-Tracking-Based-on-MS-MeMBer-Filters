function [g,chisq] = loggausspdf(xp,x0,P0)
% Calculate log of d-dimensional Gaussian prob. density evaluated at xp 
% with mean x0 and covariance P0
%
% xp: particles: dimensions d by N
% x0: mean: d by 1
% P0: covariance: d by d
%
% g: evaluation of the log-pdf
% chisq: the (xp-x)'inv(P0)(xp-x)
%
[d,N]=size(xp);
[d,n]=size(x0);
g=zeros(1,N);
%
twopi_factor = (2*pi)^(d/2);
detP_factor  = det(P0)^0.5;

if n == 1
    y=xp-x0*ones(1,N);
else
    y = xp - x0;
end;
%
chisq = sum(y.*(P0\y),1);
%
g = -(chisq/2) - log(twopi_factor) - log(detP_factor); % 1 by N

end