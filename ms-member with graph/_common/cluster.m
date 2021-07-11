function sCluster = cluster(g, vw)
% CLUSTER Decompose graph g into connected components.
% See http://blogs.mathworks.com/steve/2007/03/20/...
% connected-component-labeling-part-3/ for more.
% 函数描述：对航迹进行分簇
% 输入参数：g              航迹图
% vw            航迹图中航迹结点以及其得分值
% 输出结果：sCluster        航迹图g的连通子图
% % sCluster(iC).g : correlation matrix of each unicom component
% sCluster(iC).vw : nid in each unicom component; score of the nid
% sCluster(iC).nv : number of each unicom component
%#####输出航迹图g的连通子图sCluster
%#####vw为航迹图中航迹结点以及其得分值
% 函数说明：
% 该函数用于对航迹分簇，产生航迹图g的连通子图sCluster。
%sCluster = cluster(trackGraph, [nid; score]);
if isempty(g)
    sCluster = [];
    return
end
assert(all(diag(g)));
[p, q, r, s] = dmperm(g);%通过对这个对称矩阵进行行和列的交换，使其变成分块对角阵。则每一个分块对应图的一个connected component即簇
n = length(r)-1;
sCluster(n) = struct('g', [], 'vw', [], 'nv', []);
for iC = 1 : n
    loc = p(r(iC) : r(iC+1)-1);
    sCluster(iC).g = g(loc, loc);
    sCluster(iC).vw = vw(:, loc);
    sCluster(iC).nv = length(loc);
end

