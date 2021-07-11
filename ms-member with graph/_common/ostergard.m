function [mis, mw, sz] = ostergard(g, wl)
% OSTERGARD Use Ostergard's algorithm to find a maximum independent set.
% the boolean control variable 'found' in the original paper is problematic
%
% -- INPUT --
% g: vertex adjacency matrix
% wl: vertex weight list, all weight should be non-negative values
% 
% -- OUTPUT --
% mis: maximum independent set
% mw: total weight of mis
% sz: size of mis

    g = ~g;
    assert(~any(diag(g)));
    mw = -inf;
    nv = size(g, 1);
    c = NaN(1, nv);
    for iV = nv : -1 : 1
        clique(intersect(iV:nv, find(g(iV, :))), iV, wl(iV));
        c(iV) = mw;
        %fprintf('%d %d\n', nv+1-iV, mw);
    end

    function clique(setU, loc, w)
        if isempty(setU)
            if w > mw
                mw = w;
                mwc = loc;
                sz = length(mwc);
            end
            return
        end
        while setU
            if w + sum(wl(setU)) <= mw
                return
            end
            v = setU(1);
            if w + c(v) <= mw
                return
            end
            setU(1) = [];
            clique(intersect(setU, find(g(v, :))), [loc v], w+wl(v));
        end
    end
    mis = mwc; % max_independent_set(g) is max_weight_clique(~g)
end