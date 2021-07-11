%% slow
% function [assignments,costs]= k_best_assignment(P0,val,m)
% % get m-best ranked optimal assignment problem
% 
%   % Find the optimal and initial solution
%   [S0, C0] = mysdassignment(P0, val);
%   num_rows = size(P0,1);
%   num_cols = size(P0,2);
%   S_rows = max(P0(:,1));
% 
%   if m == 1
%     assignments(:,:,1) = S0;
%     costs = C0;
%     return; 
%   end
% 
%   % Preallocate a block of memory to hold the queue
%   N = 1000;  % Block size
%   answer_list_V = zeros(num_rows,N);
%   answer_list_S = zeros(S_rows,num_cols,N);
%   answer_list_C = NaN*ones(1,N);
% 
%   % Initialize answer list
%   answer_list_V(:,1) = val;  % cost vector for P
%   answer_list_S(:,:,1) = S0;    % solutions or assignemnts
%   answer_list_C(1) = C0;       % cost vector for problems/solutions
%   answer_index_next = 2;
% 
%   assignments = zeros(S_rows,num_cols,m);
%   costs = zeros(m,1);
% 
%   for i = 1:m
%     
%     % If all are cleared, break the loop early 
%     if all(isnan(answer_list_C))
%         max_l=answer_index_next-1;
%         assignments= assignments(:,:,1:max_l);
%         costs= costs(1:max_l,:);
%         break;
%     end
%     
%     % Grab lowest cost solution index
%     [notused,idx_top] = min(answer_list_C(1:answer_index_next-1));
%    
%     % Copy the current best solution out
%     assignments(:,:,i) = answer_list_S(:,:,idx_top);
%     costs(i,1) = answer_list_C(idx_top);
% 
%     % Copy lowest cost problem to temp
%     V_now = answer_list_V(:,idx_top);
%     S_now = answer_list_S(:,:,idx_top);
%    
%     % Delete the solution from the queue
%     answer_list_C(idx_top) = NaN;
%    
%    for a = 1:size(S_now,1)
%       % Remove it and calculate new solution
%       V_tmp = V_now;
%       [~,Loc]= ismember(S_now(a,:),P0,'rows');
%       V_tmp(Loc) = inf;
%       [S_tmp, C_tmp] = mysdassignment(P0, V_tmp);
% %     [S_tmp,C_tmp] = assignmentoptimal(V_tmp);
%     % Copy to new list
%     if size(S_tmp,1)==S_rows
%       % If we have filled the allocated space, allocate more
%       if answer_index_next > length(answer_list_C)
%         answer_list_V = cat(2,answer_list_V,zeros(num_rows,N));
%         answer_list_S = cat(3,answer_list_S,zeros(S_rows,num_cols,N));
%         answer_list_C = cat(2,answer_list_C,NaN*ones(1,N));
%       end
% 
%       answer_list_V(:,answer_index_next) = V_tmp;
%       answer_list_S(:,:,answer_index_next) = S_tmp;
%       answer_list_C(answer_index_next) = C_tmp;
%       answer_index_next = answer_index_next + 1;
% 
%       vv_tmp = V_now(Loc);
%       P_tmp = S_now(a,:);
%       P1_tmp=P0-repmat(P_tmp,size(P0,1),1);
%       idx= prod(P1_tmp,2)==0;
%       V_now(idx)=inf;
%       V_now(Loc)=vv_tmp;
%     end
%    end
%    
%   end
% 
% end

%% fast
function [assignments,costs]= k_best_assignment(P0,val,m)
% get m-best ranked optimal assignment problem

  % Find the optimal and initial solution
  [S0, C0] = mysdassignment_k_best(P0, val);
  len=min(m,length(S0));
  [~,id]=sort(C0);

  for i=1:length(len)
      assignments(:,:,i)=S0{id(i)};
      costs(i)=C0(i);
  end

end