function [sol, cost] = assignment2d(costmat)
% ASSIGNMENT Solve assignment problems involving negative cost.

func = @assignmentoptimal;
% func = @munkres;
%func = @lapjv;
% func = @assignmentoptimal;

minElement = min(costmat(:));
if minElement >= 0 % note: lapjv may be erroneous
    [sol, cost] = func(costmat);
else
    costmat = costmat - minElement;
    [sol, cost] = func(costmat);
    cost = cost + minElement*sum(sol>0);
end

% if ismember(0, sol)
%     cost = inf;
% end

if size(sol,1) > size(sol,2)
    sol = sol';
end
