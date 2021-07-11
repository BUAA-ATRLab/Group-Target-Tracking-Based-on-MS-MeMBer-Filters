function sCluster = cluster(g, vw)
% CLUSTER Decompose graph g into connected components.
% See http://blogs.mathworks.com/steve/2007/03/20/...
% connected-component-labeling-part-3/ for more.
% �����������Ժ������зִ�
% ���������g              ����ͼ
% vw            ����ͼ�к�������Լ���÷�ֵ
% ��������sCluster        ����ͼg����ͨ��ͼ
% % sCluster(iC).g : correlation matrix of each unicom component
% sCluster(iC).vw : nid in each unicom component; score of the nid
% sCluster(iC).nv : number of each unicom component
%#####�������ͼg����ͨ��ͼsCluster
%#####vwΪ����ͼ�к�������Լ���÷�ֵ
% ����˵����
% �ú������ڶԺ����ִأ���������ͼg����ͨ��ͼsCluster��
%sCluster = cluster(trackGraph, [nid; score]);
if isempty(g)
    sCluster = [];
    return
end
assert(all(diag(g)));
[p, q, r, s] = dmperm(g);%ͨ��������Գƾ�������к��еĽ�����ʹ���ɷֿ�Խ�����ÿһ���ֿ��Ӧͼ��һ��connected component����
n = length(r)-1;
sCluster(n) = struct('g', [], 'vw', [], 'nv', []);
for iC = 1 : n
    loc = p(r(iC) : r(iC+1)-1);
    sCluster(iC).g = g(loc, loc);
    sCluster(iC).vw = vw(:, loc);
    sCluster(iC).nv = length(loc);
end

