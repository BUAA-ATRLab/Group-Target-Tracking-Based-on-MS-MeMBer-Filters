function [sMwis, sBest] = graspmwis(g, vw, sGrasp)
% GRASPMWIS Find (suboptimal) maximum weighted independent sets
% See T. A. Feo et al. A Greedy Randomized Adaptive Search Procedure for
% Maximum Independent Set. Operations Reasearch.
% -- INPUT --
% g: sparse matrix showing graph adjacency 
% vw: [vertex ID list; weight value list]
% sGrasp: configure parameter structure:
%         struct('alpha', 'ntup', 'niter', 'useSearch', 'sortType')
% 
% -- OUTPUT --
% sMwis: solution structure array: struct('vertex', 'weight')
% sBest: the best solution: struct('v1', 'w1', 'v2', 'w2') for construction
%        and local search phase, respectively
%#####找出图g的独立集sMwis及最大权重独立集sBest输出给定图的独立集和最大权重独立集
% 函数描述：产生独立集，即全局假设
% 输入参数：g              航迹图
%           vw             该航迹图的航迹结点及得分
% sGrasp          求最大权重独立集问题时的参数设置
% 输出结果：sMwis          独立集，即全局假设
%           sBest           最大权重独立集，即最优假设
% 
% 函数说明：
% 找出图g的独立集和最大权重独立集，即产生全局假设和最优全局假设。
% sMwis = graspmwis(g(idx, idx), vw(:, idx), sGrasp);%#####求解多个全局假设 hypoprob(g, vw, idx, sConfig.hypoDelThresh, sConfig.sGrasp);
% assert(isequal(g, g'));
% assert(all(diag(g))); % note function rid 
% assert(all(vw(2, :) > 0));
% assert(size(g, 1) == size(vw, 2));
if issparse(g), g = full(g); end

nVertex = size(g, 1);%顶点的个数
alpha = sGrasp.alpha; % relax factor for restricted candidate list (RCL)松弛系数，取值0到1，一般为0.05
ntup = sGrasp.ntup; % maximum number of tuples
niter = sGrasp.niter; % number of randomized iteration for each tuple每一个tuple的预定的循环次数每个tuple的随机迭代次数iter
useSearch = sGrasp.useSearch; % true or false是否使用搜索算法useSearch，是为1，否则为0。
sortType = sGrasp.sortType; % 'weighted-degree' or 'weight'排序类型sortType：weight-degree，weight

% generate tuples if requested
if ntup > 0 && nVertex > 10
    nlow = min(20, nVertex); % number of vertices chosen
    sTuple = gentuple(g, vw, nlow, ntup, sortType);
    if isempty(sTuple) % no independent tuple found
        sTuple = struct('vloc', [], 'sigma', []);
    end
else
    niter = max(niter, 5);
    sTuple = struct('vloc', [], 'sigma', []);
end
ntup = numel(sTuple);

% ``GRASP'' independent sets
sMwis(ntup*niter) = struct('vertex', [], 'weight', []);
cntMwis = 0;
maxw = -inf;
for iT = 1 : ntup
    loct = sTuple(iT).vloc;
    [g1, vw1] = rid(g, vw, loct);% remove node (loc) and its neighbors from graph
    sInitSol(niter) = struct('vertex', [], 'weight', []);
    cntInitSol = 0;
    for iI = 1 : niter
        [s, w] = greedrand(g1, vw1, alpha, sortType);% find solutions using greedy randomized algorithm
        if isduplicate(s, w, cntInitSol, sInitSol)%出现重复
            continue;
        else
            cntInitSol = cntInitSol+1;
            sInitSol(cntInitSol).vertex = s;
            sInitSol(cntInitSol).weight = w;
        end
        v1 = sort([s, vw(1, loct)]);
        w1 = sum(vw(2, ismembc(vw(1,:), v1)));
        if useSearch
            [v2, w2] = localsearch(g, vw, v1);%使用局部搜索
        end
        if isduplicate(v2, w2, cntMwis, sMwis)
            continue;
        else
            cntMwis = cntMwis+1;
            sMwis(cntMwis).vertex = v2;
            sMwis(cntMwis).weight = w2;
            if w2 > maxw
                sBest = struct('v1', v1, 'w1', w1, 'v2', v2, 'w2', w2);
                maxw = w2;
            end
        end        
    end
end
sMwis = sMwis(1:cntMwis);


% ==========================================================
function sTuple = gentuple(g, vw, nlow, ntup, sortType)
% generate vertex tuples

k = 2; % tuple size
nVertex = size(g, 1);
assert(nlow <= nVertex);
sTuple = struct('vloc', {}, 'sigma', {});
if (nlow < k) || (nVertex < k)
    return
end

% sort vertices and choose nlow of them
switch sortType
  case 'weight-degree'
    [row, col] = find(g);
    % w = weighteddegree(row', col', vw(2, :));
    % [sortw, loc] = sort(w);
    w = weightratio(row', col', vw(2,:));
    [sortw, loc] = sort(w, 'descend');
    loc = loc(1:nlow);
  case 'weight'
    w = vw(2, :);
    [sortw, loc] = sort(w, 'descend');
    loc = loc(1:nlow);
  otherwise
    error('undefined sort type');
end

% generate tuples
t = 1:k;
cnt = 0;
while t
    loct = loc(t);
    if isindependent(loct, g)
        cnt = cnt + 1;
        sTuple(cnt).vloc = loct;
        %sTuple(cnt).sigma = sum(w(loct));
        [g1, w1] = rid(g, w, loct);
        sTuple(cnt).sigma = sum(w1);
    end
    t = nexttuple(t, nlow);
end

% choose best ntup tuples
sigma = [sTuple.sigma];
[sorts, loc] = sort(sigma, 'descend');
if ntup < cnt
    sTuple = sTuple(loc(1:ntup));
end


% ==========================================================
function [mis, w] = greedrand(g, vw, alpha, sortType)
% find solutions using greedy randomized algorithm

assert(alpha > 0 && alpha < 1);

mis = zeros(1, size(g, 1));
w = 0;
cnt = 0;
while ~isempty(g)
    cnt = cnt + 1;
    switch sortType
      case 'weight-degree'
        [row, col] = find(g);
        %wd = weighteddegree(row', col', vw(2, :));
        %rcl = find(wd <= (1+alpha)*min(wd));
        wr = weightratio(row', col', vw(2, :));
        rcl = find(wr >= (1-alpha)*max(wr));
      case 'weight'
        wt = vw(2, :);
        rcl = find(wt >= (1-alpha)*max(wt));
      otherwise
        error('undefined sort type');
    end
    
    % randomly choose one element of rcl
    loc = rcl(unidrnd(numel(rcl)));%产生从1到numel(rcl)所指定的最大数数之间的离散均匀随机整数
    mis(cnt) = vw(1, loc);
    w = w + vw(2, loc);
    
    [g, vw] = rid(g, vw, loc);
end
mis = sort(mis(1:cnt));


% ==========================================================
function [sol, wsol] = localsearch(g, vw, s)
% perform local search to improve working solution s

v = vw(1, :);
w = vw(2, :);
sol = s;
%wsol = sum(w(ismember(v, s)));
assert(issorted(s));
wsol = sum(w(ismembc(v, s)));
ns = length(s); % number of vertices in the solution 
nv = length(v); % number of vertices in the graph

% search for local optimization
for k = 1 : ns
    idx = true(1,ns); 
    idx(k) = false;
    ss = s(idx);
    
    %tf = ismember(v, ss);
    tf = ismembc(v, ss);
    wss = sum(w(tf));
    loc = find(tf);
    [g1, vw1] = rid(g, vw, loc);
    
    if ~isempty(g1)
        [loc, wm, sz] = ostergard(g1, vw1(2, :));
        mis = vw1(1, loc(1:sz));
        if wm + wss > wsol
            sol = sort([ss mis]);
            wsol = wm + wss;
        end
    end
end


% ==========================================================
function tuple = nexttuple(tuple, n)
% return a subset of 1:n next to input tuple in lexical order

nt = length(tuple);
loc = nt;
while loc && tuple(loc) == n-(nt-loc)
    loc = loc - 1;
end
if loc == 0
    tuple = [];
else
    tuple(loc:end) = tuple(loc) + (1:(nt-loc+1));
end


% ==========================================================
function tf = isindependent(loc, g)
% check for vertex independency
gg = tril(g(loc, loc), -1);
tf = ~any(gg(:));


% ==========================================================
function tf = isduplicate(s, w, cnt, sMwis)
% true for duplicate solution重复的
tf = false;
if cnt == 0, return; end
tol = 1e-6; % tolerance
weight = [sMwis.weight];
locCand = find(abs(weight-w) < tol);
if isempty(locCand), return; end
for iC = 1 : length(locCand)
    cand = sMwis(locCand(iC)).vertex;
    assert(issorted(s));
    assert(issorted(cand));
    if isequal(s, cand)
        tf = true;
        return
    end
end


% ==========================================================
function [g, vw] = rid(g, vw, loc)
% remove node (loc) and its neighbors from graph

if isempty(loc)
    return
elseif length(loc) == 1
    %rm = [loc find(g(loc, :))];
    rm = false(1, size(g, 1));
    rm(loc) = true;
    rm = or(rm, g(loc, :));
else
    %s = sum(g(loc, :));
    %rm = [loc find(s)];
    rm = false(1, size(g, 1));
    rm(loc) = true;
    for iL = 1 : length(loc)
        rm = or(rm, g(loc(iL), :));
    end
end

kp = ~rm;
g = g(kp, kp);
vw = vw(:, kp);