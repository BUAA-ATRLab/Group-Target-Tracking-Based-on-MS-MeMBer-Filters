% auction.m 
% 2-D auction algorithm, modified version of Yanhua Ruan's auc.m
% see also: auction_2d.m 
%           for 2-D auction dealing with dummy measurements (missed detections)
function [q,omiga,assign] = auction(cost)

assert(max(cost(:)) < inf);
valid = isfinite(cost);
v = cost(valid);
small = -sum(abs(v))-1;
cost(~valid) = small;

[m, n] = size(cost);
assign = zeros(m,1);

if m > n
    [q,omiga,assign_t]=auc(cost');
    omiga = omiga';

    for ind=1:m
        if assign_t(ind) > 0
            assign(assign_t(ind)) = ind;
        end
    end
else
    [q,omiga,assign]=auc(cost);
end

if q < sum(v(v < 0))
    q = -inf;
end
assign = assign';



% Auction algorithm for solving 2-D Assignment problem
% by Yanhua Ruan
% June.12,1999, modified in Dec.20,1999, a Gauss-Seidel version
function [q,omiga,assign]=auc(befi)

% AUCTION_2D_SCALE=0.25;
% AUCTION_2D_MAXC=1000.0;
% AUCTION_2D_FACTOR=8;
% AUCTION_2D_DECRFACTOR=4;  % 4<=. . <=10
% AUCTION_2D_ENDEPS=1;
AUCTION_2D_MAXTHRESH=100;
% AUCTION_2D_SMALL=realmin;


% scaling, make 1. all elements be positive; 2. the smallest element be 1.
s1=min(befi(:));
if s1<0
    befi=befi-s1;
end
[nu1,nu2,v]=find(befi);
s2=min(v); % spz befi>0 already. 
befi=befi./s2;

%  init. price for each object
[n_person,n_object]=size(befi);
price=zeros(1,n_object); 

assign=zeros(n_object,1);
curlst=linspace(1,n_person,n_person);
numnew = n_person;

thresh = min(n_person/5,AUCTION_2D_MAXTHRESH);
epsilon = 1e-6;
decr=epsilon;

while (numnew>thresh)||((numnew>0))
    nolist=numnew;  % the number of unsigned person
    numnew=0;

    for i=1:nolist
        row = curlst(i);
        
		tmax=befi(row,1:n_object)-price(1:n_object);
		[max1 bstcol]=max(tmax);
		tmax(bstcol)=-realmax;
		max2=max(tmax);
        
        if bstcol ~= 0
            price(bstcol) = price(bstcol)+max1-max2+decr; 
            oldrow = assign(bstcol);
            assign(bstcol) = row;
            if oldrow > 0
                numnew = numnew + 1;            
                curlst(numnew) = oldrow;
            end
        end   
    end 
    %decr=max(decr/AUCTION_2D_DECRFACTOR,epsilon);
    %decr=max(decr/1.414, .0001);
end  % while

omiga=zeros(n_person,n_object);

q=0;
for i=1:n_object
    if assign(i)~=0
        omiga(assign(i),i)=1;
        q=q+befi(assign(i),i);      
    end   
end   
% scaling back
q = q*s2+s1*n_person*(s1<0);

