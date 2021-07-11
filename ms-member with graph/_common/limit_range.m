function clipped_r = limit_range(r)
% clips the value of r to within [lower, upper]. Usefull for clipping the
% values of probabilities for numerical stability. 
% Input:  r - array or cell or array of probabilities.
% Output: clipped_r - same as r but with values clipped at upper and
%         lower.

upper = 1 - 1e-3;
lower = 1e-10;

% Clip the existence probabilities within [lower, upper]
if ~iscell(r)   % if r is an array of double
    r(r > upper) = upper;
    r(r < lower) = lower;
    
    clipped_r = r;
else            % if r is an array of cells
    [r{[r{:}] > upper}] = deal(upper);
    [r{[r{:}] < lower}] = deal(lower);
    
    clipped_r = r;
    
end


end