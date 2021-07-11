function[dik]= mahalanobis_distance(Xi,Xk,Pi,Pk)
%mahalanobis_distance.m
% This program calculates the mahalanobis distance between two vectors
%
% Inputs : 
%         Xi, Xk are the target's states e.g Xi=[x,y,vx,vy]
%         Pi, Pk are the target's variance-covariance 3D matrices 


% Output :
%         dik is the mahalanobis distance between Xi and Xk
%----Amadou Dec 2007----

n=length(Xi);
Pik=zeros(n,n);

for j=1:n
    Pik(j,j)=Pi(j,j)+Pk(j,j);
end

dik=sqrt((Xi-Xk)*inv (Pik)*(Xi-Xk)');