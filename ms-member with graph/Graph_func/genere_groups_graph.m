function[Adj]= genere_groups_graph(X,P,epsilonx,epsilonv)
%genere_groups_graph.m
%  This program Generates the graph associated to the groups of targets
%
% Inputs : 
%         X is the target's states matrix where X(:,i)=[x,y,vx,vy]' is
%         vector of the state.  
%         P is the target's variance-covariance 3D matrix where P(:,:,i) is
%         a 4*4 matrix
%         epsilon is a treshold to decide weither two targets are in the
%         same group or not,
% Output :
%         Adj is the adjacent matrix associated to the edges of the
%         undirected graph of the group structure.  

%----Amadou Dec 2007----


[~,N]=size(X);
Adj =zeros(N,N);


for i=1:N-1
    for k= i+1:N
        dik1=mahalanobis_distance([X(1,i) X(2,i)],[X(1,k) X(2,k)],P([1 2],[1 2],i),P([1 2],[1 2],k));
        dik2=mahalanobis_distance([X(3,i) X(4,i)],[X(3,k) X(4,k)],P([3 4],[3 4],i),P([3 4],[3 4],k));
        if (dik1 <= epsilonx) && (dik2 <= epsilonv)%(dik2 <= 5)
            Adj(i,k)=1;
            Adj(k,i)=1;% added in order to use the "torsche" toolbox(graph theory)
        end
    end
end