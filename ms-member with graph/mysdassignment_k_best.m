

function [gam_all, cost] = mysdassignment_k_best(sub, val)

USE_AUC = 0;
gam_all=cell(1,0);
n = max(sub);
S = size(sub, 2);
assert(S >= 2);
if S==2
    gam_all{1}=[1:n(1)]';
    c = buildmatrix(sub, val, n(1), n(2), 'full');
    if USE_AUC
        [J,~,assign]=auction(-c);
        cost(1)=-J;
    else
        [assign,J]=assignment2d(c);
        cost(1)=J;
    end
    assi=find_assi_index(assign, n(1));
    gam_all{1}=[gam_all{1} assign'];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_dual=-realmax;
f_primal=realmax;
rgap=realmax;
iter=0;
mingap=.01;  %  typically 0.01 to 0.05
maxiter=500;
igam=0;

a=.1;   % typically, .05<= a <=.3
b=1.2;  % typically, 1.1<= b <=1.6
beta=1;
alpha=2;
g=ones(S,max(n));
u=zeros(S,max(n));
H=zeros(max(n),max(n),S);

while rgap>mingap && iter<maxiter
    gam=[1:n(1)]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2: Compute reduced costs for r=S-1, S-2,...,2.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2A: c-u
    for r=2:S-1 
        val2 = val;
        for i = 1 : size(sub, 1)
            ith = sub(i, :);
            ind = sub2ind(size(u), r+1:S, ith(r+1:S));
            val2(i) = val(i) - sum(u(ind));
        end
        
        nr = prod(n(1:r));
        nc = prod(n)/nr;
        c2 = buildmatrix(sub, val2, nr, nc, 'mysparse');
        sub2 = zeros(numel(c2), 2);
        val2 = zeros(numel(c2), 1);
        for i = 1 : numel(c2)
            m = c2{i};
            [x, j] = min(m(:, 3));
            sub2(i, :) = [m(j,1) m(j,2)];
            val2(i) = x;
        end
        di = sub2;
        clear c2;
        %disp(di);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STEP 3: Enforce constraint set r. r=2,3,...,S
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % !!!! some work need to do for index (index(r))?
        [d2, ir] = buildmatrix(sub2(:,1), val2, prod(n(1:r-1)), n(r), 'mysparse');
        uir = unique(ir);
        dr = zeros(n(1), n(r));
        for i=1:n(1)
            index=gam(i,r-1)-1;		 % gam(i,1) = i
            for j=r-2:-1:1
                index=index*n(j)+gam(i,j)-1;
            end   
            index=index+1;
            loc = find(index == uir);
            if ~isempty(loc)
                m = d2{loc};
                dr(i, m(:,2)) = m(:,3);
            end
        end
        dr(dr == 0) = inf;
        
        if USE_AUC
            %  dr is always a two-dim (n1 X n(r)) matrix
            [J(r),omiga,assign]=auction(-dr); 
            J(r)=-J(r)+sum(sum(u(r+1:S,:)));
        else
            [assign,J(r)]=assignment2d(dr);
            J(r)=J(r)+sum(sum(u(r+1:S,:)));
        end
        % disp(J(r));
        % disp(assign);
        % disp('================');
        clear d2 dr;
        assi = find_assi_index(assign, n(1));
        gam=[gam assign'];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STEP 4: Update Lagrangian Multiplers u for r=3,4,...,S.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if r==2
            Jstar=J(2);   
            f_dual=max(f_dual,Jstar);%报告中式（1.15）
        else    
            g(r,1:n(r))=ones(1,n(r));
            for k=1:n(1)
                index=gam(k,r-1)-1;		 % gam(k,1) = k
                for i=r-2:-1:1
                    index=index*n(i)+gam(k,i)-1;
                end
                index=index+1;
                loc = find(index == do(:,1));
                if ~isempty(loc)
                    j=mod(do(loc,2)-1,n(r))+1; % ????
                    j = round(j);
                    g(r,j)=g(r,j)-1;
                end
            end
            
            % find the adapative factor for the updating of u
            % Accelerated Subgradient Update:
            % adapt=eye(n(r));
            % u(r,1:n(r))=u(r,1:n(r))+g(r,1:n(r))*adapt;
            if g(r,1:n(r))==0
            else
                if mod(iter,n(r))==0
                    p=zeros(n(r),1);
                    H2=eye(n(r));
                else
                    H2=H(1:n(r),1:n(r),r);
                    p=H2*g(r,1:n(r))';
                    H2=H2+(1-alpha^(-2))*p*p'/(g(r,1:n(r))*p);%报告中式（1.20）
                end   
                fa_dual=(1+a/beta^b)*f_dual;
                adapt=(alpha+1)/alpha*(J(r)-fa_dual)/norm(g(r,1:n(r)))/norm(g(r,1:n(r)));
                u(r,1:n(r))=u(r,1:n(r))+p'*adapt;%报告中式（1.19）
                if J(2)<f_dual
                    beta=beta+1;%报告中式（1.21）
                else
                    beta=max(beta-1,1);
                end   
                H(1:n(r),1:n(r),r)=H2;
            end 
        end
        do = di;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STEP 5: Recursion: Identify assignments for the R-D problem
        %         Increment R and go to the STEP 3.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end % r or R
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 6: Iteration: Improve solution quality, go to STEP 2 or STEP 7.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
    % for r=S case.
    %d2=reshape(c,[prod(n(1:S-1)),n(S)]);
    [d2, ir] = buildmatrix(sub, val, prod(n(1:S-1)), n(S), 'mysparse');
    uir = unique(ir);
    dr = zeros(n(1), n(S));
    for i=1:n(1)
        index=gam(i,S-1)-1;		 % gam(i,1) = i
        for j=S-2:-1:1
            index=index*n(j)+gam(i,j)-1;
        end
        index=index+1;
        loc = find(index == uir);
        if ~isempty(loc)
            m = d2{loc};
            dr(i, m(:,2)) = m(:,3);
        end
    end
    dr(dr == 0) = inf;
    %save('mysda.mat', 'c', 'd2', 'dr', 'gam'); break;
    
    if USE_AUC
        [J(S),omiga,assign]=auction(-dr);
        J(S)=-J(S);
    else
        [assign,J(S)]=assignment2d(dr);
    end
    % disp(J(S));
    % disp(assign);
    % disp('================');
    clear d2 dr;

    assi=find_assi_index(assign,n(1));
    gam=[gam assign'];
    
    % update multiplier u(S,:) 
    g(S,1:n(S))=ones(1,n(S));
    for k=1:n(1)
        index=gam(k,S-1)-1;		 % gam(k,1) = k
        for i=S-2:-1:1
            index=index*n(i)+gam(k,i)-1;
        end
        index=index+1;
        loc = find(index == do(:,1));
        if ~isempty(loc)
            j=mod(do(loc,2)-1,n(S))+1; % ????
            j = round(j);
            g(S,j)=g(S,j)-1;
        end
    end
    
    if g(S,1:n(S))==0
    else   
        if mod(iter,n(S))==0
            p=zeros(n(S),1);
            H2=eye(n(S));
        else  
            H2=H(1:n(S),1:n(S),S);
            p=H2*g(S,1:n(S))';
            H2=H2+(1-alpha^(-2))*p*p'/(g(S,1:n(S))*p);%报告中式（1.20）
        end   
        fa_dual=(1+a/beta^b)*f_dual;
        adapt=(alpha+1)/alpha*(J(S)-fa_dual)/norm(g(S,1:n(S)))/norm(g(S,1:n(S)));
        u(S,1:n(S))=u(S,1:n(S))+p'*adapt;      %报告中式（1.19）  
        if J(2)<f_dual
            beta=beta+1;%报告中式（1.21）
        else
            beta=max(beta-1,1);
        end   
        H(1:n(S),1:n(S),S)=H2;
    end   
    
    % Jstar=J(2);   
    % f_dual=max(f_dual,Jstar);
    f_primal=min(f_primal,J(S));%报告中式（1.16）
    rgap=(f_primal-f_dual)/abs(f_primal);%报告中式（1.17）
    iter=iter+1;%报告中式（1.18）
    
    %update gam
    flag=0;
    for ii=1:igam
        if isequal(gam_all{ii},gam)
            flag=1;
            break;
        end
    end
    if flag==0
        igam=igam+1;
        gam_all{igam}=gam;
    end
    
end  % while

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 7: Final Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the last assignment result, make it easy to find obj(s(i)) for person(i)

%assignlast=gam(:,2);
%assignlast2=gam(:,3);
%[temp, assignlast]=sort(assign);
%[temp, assignlast2]=sort(assign2);
%assignlast=assignlast(1+n2_FA:n2);  
%assignlast2=assignlast2(1+n3_FA:n3);
for i=1:igam
    v = all(gam_all{i}, 2);
    gam_all{i} = gam_all{i}(v, :);
    c1 = mat2cell(sub, size(sub,1), ones(1,size(sub,2)));
    ind1 = sub2ind(n, c1{:});
    c2 = mat2cell(gam_all{i}, size(gam_all{i},1), ones(1, S));
    ind2 = sub2ind(n, c2{:});
    [tf, loc] = ismember(ind2, ind1);
    gam_all{i} = gam_all{i}(tf, :);
    cost(i) = sum(val(loc(tf)));
end






% function assi = find_assi_index(assign, n)
%
% assign  ----- assignment vector indicating which measurement goes to which track
%               assign(j) = k means the j-th measurement goes to track k
%               assign(j) = 0 means the j-th measurement is unassigned
% n       ----- number of tracks
% assi    ----- assignment index
%               assin(k) = j means the k-th track associates j-th measurement
function assi = find_assi_index(assign, n)

assi = zeros(n, 1);
for i=1:n
    tmp = find(~(assign-i));
    if ~isempty(tmp)
        assi(i) = tmp;
    end
end



% build 2d sparse/full matrix according to input subscripts and values
function [m, ind1, ind2] = buildmatrix(sub, val, nRow, nCol, type)

assert(size(sub, 1) == length(val));

n = size(sub, 1);
s = size(sub, 2);
if s == 1
    sub = [sub ones(n,1)];
    s = 2;
end
dim = max(sub);
c = mat2cell(sub, n, ones(1, s));
ind = sub2ind(dim, c{:});
switch type
  case 'full'
    m = zeros(nRow, nCol);
    m(ind) = val;
  case 'sparse'
    [ind1, ind2] = ind2sub([nRow, nCol], ind);
    m = sparse(ind1, ind2, val, nRow, nCol);
  case 'mysparse'
%      ind1 = sub2ind(dim(1:end-1), c{1:end-1});
    [ind1, ind2] = ind2sub([nRow, nCol], ind);
    ijs = sortrows([ind1(:), ind2(:), val(:)], 1);
    [uids, uidx] = unique(ijs(:,1), 'last');
    sublen = [uidx(1); diff(uidx)];
    m = mat2cell(ijs, sublen, 3);
  otherwise
    error(['unknown output type ''' type '''']);
end