function[G,Adj_O]= genere_groups_centre_graph(X,P,G,epsilon)
%genere_groups_centre_graph.m
% This program Generates the graph associated to the centres of groups
%
% Inputs : 
%         X is the target's states matrix where X(i,:)=[x,y,vx,vy] is
%         vector of the state,  
%         P is the target's variance-covariance 3D matrix where P(i,:,:) is
%         a 4*4 matrix,
%         epsilon is a treshold to decide weither two groups are
%         neighbouring or not,
%         Adj is the edges of the graph constitued by the targets.

% Output :
%         G is a structure, with the group centre, the indexes of targets i
%         in the groups and the variance associated to each centre,
%         Adj_o is the  adjacent matrix associated to the edges of the
%         undirected graph constituted by the groups centres.  

%----Amadou Dec 2007----
%G=connected_components(Adj);
nG=length(G);
n=length(X(1,:));%4
for i=1:nG
    G(i).centre =mean(X(G(i).index,:),1);
    ng=length(G(i).index);
    G(i).P= zeros(1,n,n);
    for k=1:ng
        G(i).P=G(i).P+P(G(i).index(k),:,:);
    end
    G(i).P=G(i).P/ng;
end
%Edge creation
Adj_O =zeros(nG,nG);
for i=1:nG-1
    for k= i+1:nG
        dik=mahalanobis_distance(G(i).centre,G(k).centre,G(i).P,G(k).P);
        if dik <= epsilon
            Adj_O(i,k)=1;
            Adj_O(k,i)=1;
        end
    end
end