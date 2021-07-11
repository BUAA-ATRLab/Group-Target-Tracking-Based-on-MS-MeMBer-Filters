function wd = weightratio(row, col, wl)
% WEIGHTRATIO Calculate an array of W(v)/sum(weight of N(v))

wd = zeros(1, length(wl));
for loc = 1 : length(row)
    vertex = col(loc);
    neighbor = row(loc);
    wd(vertex) = wd(vertex) + wl(neighbor);
end
wd = wl./wd;

