function logsum = logsumexp(w, varargin)

%performs log-sum-exp trick to avoid numerical underflow
% input:  w weight vector or matrix assumed already log transformed
%         dim (optional) for matrices w specifies the dimension on which to
%                        perform logsumexp
% output: log(sum(exp(w)))

%[nb_rows, nb_cols] = size(w); 

%if nb_rows*nb_cols == length(w)    % w is a vector

    if all(w ==-inf)
        logsum= -inf;
        return;
    end
    
    [val,~] = max(w);
    logsum = log(sum(exp(w-val))) + val;
    
% else  % w is a matrix, perform logsumexp along dimension dim 
%     
%     if ~isempty(varargin)
%         dim = varargin{1};
%     end
%     
%    [val,~] = max(w, [], dim); 
%    
%    if dim == 1, 
%        valm = repmat(val,nb_rows,1);
%    else
%        valm = repmat(val,1, nb_cols);
%    end
%    
%    logsum = log(sum(exp(w-valm), dim)) + val;
%    
%    logsum( all(w == -Inf, dim) ) = -Inf;
%    
% end

end
