function[G]= connected_components(E_matrix)
%connected_components.m
% This program Generates the groupe structure G using the adjacent matrix
%
% Inputs : 
%         E_matrix is the graph adgacent matrix

% Output :
%        G is the set of groups
% code from
% http://blogs.mathworks.com/steve/2007/03/20/connected-component-labeling-part-3/
if isempty(E_matrix)
   G=[];
   return;
end
n= length(E_matrix);
E_matrix(1:n+1:end) = 1;
[p,q,r,s] = dmperm(E_matrix);
nG=length(r);
for i=1:nG-1
    %G(i).index=p(r(i):r(i+1)-1);
    G(i).index=p(r(i+1)-1:-1:r(i));
end 


% 
% g=graph('adj',E_matrix);
% CC=tarjan(g);
% nG=max(CC);
% 
% for i=1:nG
%     G(i).index=find(CC==i);
% end   