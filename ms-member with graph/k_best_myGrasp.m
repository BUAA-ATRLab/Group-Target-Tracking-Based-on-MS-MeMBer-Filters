%% 最大权重独立集求解
function [assignments,costs]= k_best_myGrasp(P0,val,m,sGrasp)
% get m-best solutions

  % Find the optimal and initial solution
row = size(P0,1);
mat=zeros(row,row);
for i=1:row
    temp=repmat(P0(i,:),row,1)-P0;
    idx= prod(temp,2)==0;
    mat(i,idx)=1;
end
trackGraph=sparse(mat);
sCluster = cluster(trackGraph, [1:row;val(:)']);


for iC=1:length(sCluster)
sMwis{iC} = graspmwis(full(sCluster(iC).g), sCluster(iC).vw, sGrasp);
end
id_mat=cell2mat({sMwis{1}.vertex}');
w_mat=[sMwis{1}.weight]';
min_len=min(m,size(id_mat,1));
[w_mat, loc] = sort(w_mat, 'descend');
w_mat=w_mat(1:min_len);
id_mat=id_mat(loc(1:min_len),:);
for i=2:length(sMwis)
    cHypo = cell2mat({sMwis{i}.vertex}');%#####全局假设
    hypoScore = [sMwis{i}.weight]';%#####相应的假设得分
    temp_S=kron(hypoScore,ones(size(id_mat,1),1));
    temp_h=kron(cHypo,ones(size(id_mat,1),1));
    id_mat=[repmat(id_mat,size(cHypo,1),1) temp_h];
    w_mat=repmat(w_mat,size(cHypo,1),1)+temp_S;
    min_len=min(m,size(id_mat,1));
   [w_mat, loc] = sort(w_mat, 'descend');
    w_mat=w_mat(1:min_len);
   id_mat=id_mat(loc(1:min_len),:);
end
min_len=min(m,size(id_mat,1));
[w_mat, loc] = sort(w_mat, 'descend');
w_mat=w_mat(1:min_len);
id_mat=id_mat(loc(1:min_len),:);
for i=1:min_len
assignments(:,:,i)=sortrows(P0(id_mat(i,:),:));

costs(i)=sum(val(id_mat(i,:)));
end
end