function output_member = filter_ms_member_kf(filter, model, meas, truth)


%output variables
output_member.state_X_hat = cell(meas.K,1);   % inferred target states
output_member.state_n_hat = zeros(meas.K,1);  % estimated number of targets
output_member.filter      = filter;           % store filtering params

%=== Filtering
member_mark_start = tic;      % start timing the filter

%Graph
Graph=cell(meas.K,1);

% void initial prior
MS_MeMBer_update.M = 0;       % nb of Bernoulli components
MS_MeMBer_update.r = [];      % prob of existence of the Bernoulli components
MS_MeMBer_update.m = [];      % means of Bernoulli components
MS_MeMBer_update.P = [];      % covariance matrices of Bernoulli components
MS_MeMBer_update.G = [];      % nb of graph structure(0 is single target without group)

%recursive filtering
for k = 1:meas.K
    
    % Prediction and update
    MS_MeMBer_predict = predict(MS_MeMBer_update, model); 
    
    MS_MeMBer_update  = update(MS_MeMBer_predict, model, filter, meas.Z{k});                      
    
    % Pruning 
    MS_MeMBer_update = prune(MS_MeMBer_update, filter);                                            
    
    % Cap the Bernoulli components at either the maximum number of allowed
    % components or to a number proportional with the number of estimated
    % targets
    N_cap = filter.comp_max;
    MS_MeMBer_update = cap(MS_MeMBer_update, N_cap);                                             
    % State estimation
    [output_member.state_X_hat{k}, output_member.state_n_hat(k), ind] = extract_estimates(MS_MeMBer_update);
    % Return the updated Bernoulli existence probabilities for diagnostics
    output_member.rkk{k} = MS_MeMBer_update.r  / sum(MS_MeMBer_update.r);   

end

% return filter computation time 
output_member.all_time    = toc(member_mark_start);

% calculate and display average OSPA error
error_t = zeros(1,truth.K);
for k = 1:truth.K
    error_t(k) = ospa_dist( get_comps(truth.X{k},[1 2]), get_comps(output_member.state_X_hat{k},[1 2]), model.ospa_cutoff, model.ospa_order);
end

% return OSPA errors
output_member.ospa = error_t;
% return structure
output_member.Graph = Graph;  
end % end function display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function MS_MeMBer_predict = predict(MS_MeMBer_update, model)
% Predict PHD function
% Predict surviving tracks
MS_MeMBer_predict   = MS_MeMBer_update;
MS_MeMBer_predict.r = MS_MeMBer_update.r * model.kin.pS;
if MS_MeMBer_update.M
    [MS_MeMBer_predict.m, MS_MeMBer_predict.P] = kalman_predict_multiple(model, MS_MeMBer_update.m, MS_MeMBer_update.P);
end


%% Sort the surviving components in oder of decreasing weights. In the case of overlapping subsets for two predicted Bernoulli components, the partition formation procedure assigns a measurement subset to the first predicted component and a miss detection for the second. Hence, it is better to have the highest component first.  
if MS_MeMBer_predict.M
   [~, sorted] = sort(MS_MeMBer_predict.r,'descend');

   MS_MeMBer_predict.r = MS_MeMBer_predict.r(sorted);
    
    for i = 1:MS_MeMBer_predict.M
        sorted_m(:,i) = MS_MeMBer_predict.m(:,sorted(i));
        sorted_P(:,:,i) = MS_MeMBer_predict.P(:,:,sorted(i));
    end
    MS_MeMBer_predict.m = sorted_m;
    MS_MeMBer_predict.P = sorted_P;  
end

% append birth components
for i = 1:model.birth.J_gamma
    MS_MeMBer_predict.r(MS_MeMBer_predict.M + i) = model.birth.r_gamma(i);
    MS_MeMBer_predict.m(:,MS_MeMBer_predict.M + i) = model.birth.m_gamma(:,i);
    MS_MeMBer_predict.P(:,:,MS_MeMBer_predict.M + i) = model.birth.P_gamma;
end

% account for the increase in the components
MS_MeMBer_predict.M = MS_MeMBer_update.M + model.birth.J_gamma;

end

function MS_MeMBer_predict = predict1(MS_MeMBer_update, model, Gr)
% Predict PHD function
% Predict surviving tracks
MS_MeMBer_predict   = MS_MeMBer_update;
MS_MeMBer_predict.r = MS_MeMBer_update.r * model.kin.pS;
if MS_MeMBer_update.M
    [MS_MeMBer_predict.m, MS_MeMBer_predict.P] = kalman_predict_multiple1(model, MS_MeMBer_update.m, MS_MeMBer_update.P, MS_MeMBer_update.G, Gr);
end


%% Sort the surviving components in oder of decreasing weights. In the case of overlapping subsets for two predicted Bernoulli components, the partition formation procedure assigns a measurement subset to the first predicted component and a miss detection for the second. Hence, it is better to have the highest component first.  
if MS_MeMBer_predict.M
   [~, sorted] = sort(MS_MeMBer_predict.r,'descend');

   MS_MeMBer_predict.r = MS_MeMBer_predict.r(sorted);
    
    for i = 1:MS_MeMBer_predict.M
        sorted_m(:,i) = MS_MeMBer_predict.m(:,sorted(i));
        sorted_P(:,:,i) = MS_MeMBer_predict.P(:,:,sorted(i));
    end
    MS_MeMBer_predict.m = sorted_m;
    MS_MeMBer_predict.P = sorted_P;  
end

% append birth components
for i = 1:model.birth.J_gamma
    MS_MeMBer_predict.r(MS_MeMBer_predict.M + i) = model.birth.r_gamma(i);
    MS_MeMBer_predict.m(:,MS_MeMBer_predict.M + i) = model.birth.m_gamma(:,i);
    MS_MeMBer_predict.P(:,:,MS_MeMBer_predict.M + i) = model.birth.P_gamma;
end

% account for the increase in the components
MS_MeMBer_predict.M = MS_MeMBer_update.M + model.birth.J_gamma;

end

function MS_MeMBer_update = update(MS_MeMBer_predict, model, filter, Z)
% Update the predicted PHD function and cardinality distribution
% Check to see if the current measurement set is not empty
empty_data = 1;
for s = 1:model.obs.nb_sensors
    if ~isempty(Z{s})
        empty_data = 0;
        break;
    end
end


gamma = prod(1-model.obs.pD);  % probability that all sensors miss the targets
    
if empty_data
    % if no sensor measurements are available, we keep the predicted
    % Bernoulli components but with decreased wights.
 
    MS_MeMBer_update = MS_MeMBer_predict;
    MS_MeMBer_update.r = gamma * MS_MeMBer_predict.r ./ (1 - MS_MeMBer_predict.r +  MS_MeMBer_predict.r * gamma);
else      
    %% Find best observation subsets (paths) for each PHD component
   [trellis_paths, m_paths, P_paths, dW] = find_best_subsets(model, filter, MS_MeMBer_predict, Z);
  
    %% Find likely partitions
   [~, unique_partitions, alpha_P] = find_best_partitions(MS_MeMBer_predict.M, trellis_paths, dW, filter.P_max, model.obs.lambda_c);
  
    index = 0;                                                             % counter for updated Bernoulli components 
    r_kk = [];
    
    for j = 1:MS_MeMBer_predict.M                                          % for each predicted Bernoulli component 
        
        subsets_j = unique_partitions(:,j);                                % select subsets assigned to the j-th Bernoulli component     
        [unique_subsets_j] = unique(subsets_j,'rows');                     % consider only the unique subsets
        
        nb_unique_subsets_j = length(unique_subsets_j);
               
        for  i = 1:nb_unique_subsets_j                                     % an updated Bernoulli component is created for each unique subset 
            index = index + 1; 
            
            if unique_subsets_j(i) ~= 0                                     % if subset is not the empty subset \emptyset_{1:s}
                r_kk(index) = sum( alpha_P( subsets_j ==  unique_subsets_j(i) ) );  % add the weight of all partitions that assign the subset unique_subsets_j(i) to the j-th Bernoulli
                MS_MeMBer_update.r(index) = limit_range(r_kk(index));               % clip the component prob. of existane to 0.999 (for numerical stability when calculating cardianlity pmf) 
                                
                MS_MeMBer_update.m(:,index) = m_paths{j}(:,unique_subsets_j(i));    % updated mean and covariance                 
                MS_MeMBer_update.P(:,:,index) = P_paths{j}(:,:,unique_subsets_j(i));
            else                                                                     % if subset is the empty subset => legacy component, use the predicted Bernoulli and only update the existece probability
                r_kk(index) = sum( alpha_P( subsets_j ==  0 ) ) * gamma * MS_MeMBer_predict.r(j) ./ (1 - MS_MeMBer_predict.r(j) +  MS_MeMBer_predict.r(j) * gamma);
                MS_MeMBer_update.r(index) = limit_range(r_kk(index));
                
                MS_MeMBer_update.m(:,index) = MS_MeMBer_predict.m(:,j);    % mean and covaraicne are unchanged               
                MS_MeMBer_update.P(:,:,index) = MS_MeMBer_predict.P(:,:,j);
             end
        end       
    end
    
    MS_MeMBer_update.M = length(MS_MeMBer_update.r);

end  % if ~empty data
end % end function update


function [best_paths, m_paths, P_paths, dW] = find_best_subsets(model, filter, MS_MeMBer_predict, Z)    
% Function that searches for the best W_max subsets for each PHD component

M_kk1 = MS_MeMBer_predict.M; % number of predicted PHD components
gamma = prod(1 - model.obs.pD);
global flag
% initialize output variables
best_paths = cell(1,M_kk1);
dW = cell(1,M_kk1);

% for each component
for j = 1:M_kk1

    curr_J  = 1; % the number of current subsets W, starts with 1 empty W
    curr_paths = []; % the paths of all the current subsets W
    % we initilize the path parameters with the predicted parameters 
    curr_w  = MS_MeMBer_predict.r(j) ; % we assume the mass of the current component as initial subset weight 
    curr_m  = MS_MeMBer_predict.m(:,j);
    curr_P  = MS_MeMBer_predict.P(:,:,j);
  
   % process each sensor sequentially
    for s = 1:model.obs.nb_sensors
        curr_Z = Z{s};
        nb_meas = size(curr_Z, 2);    % current number of measurements   

        % Kalman updates for the for no detected observation case
        local_w = curr_w*(1 - model.obs.pD(s));
        local_m = curr_m;
        local_P = curr_P;
        
        if nb_meas ~= 0   % if there is at least one measurement
          
            [qz_update, m_update, P_update] = kalman_update_multiple(curr_Z, model, curr_m, curr_P);
           
            if flag==1
                new_w = (model.obs.pD(s)/model.obs.pdf_c/model.obs.lambda_c*repmat(curr_w, nb_meas,1) .* qz_update')';
            else
                new_w = (model.obs.pD(s)/model.obs.pdf_c*repmat(curr_w, nb_meas,1) .* qz_update')';
            end
            
            local_w = [local_w new_w(:)'];
            local_m = [local_m, m_update];
            local_P = cat(3, local_P, repmat(P_update, 1,1, nb_meas));
        end
 
        local_paths = repmat(curr_paths, nb_meas+1, 1);
        new_obs     = kron(1:nb_meas, ones(curr_J, 1));
        new_obs     = [-1*ones(curr_J,1); new_obs(:)];
        local_paths = [local_paths new_obs];
        

        % select best observation sequences till now 
        % always include the empty observation path first
        [~,index] = sort(local_w(2:end),'descend');
        % select only the best W_max measurements as candidates for the j component
        J = min(length(index), filter.W_max);  % update number of components propagated
        index = index(1:J) + 1; % account for the void measurement association
        curr_w = [local_w(1) local_w(index)]; % first association is with the void measurement
        curr_m = [local_m(:,1) local_m(:,index)];
        curr_P = cat(3, local_P(:,:,1), local_P(:,:,index));
        curr_paths = [local_paths(1,:); local_paths(index,:)];
        curr_J = J + 1;
    end % for each sensor
    
    % final paths information for j-th Gaussian component
    best_paths{j} = curr_paths;
    m_paths{j}    = curr_m;
    P_paths{j}    = curr_P;
    
    curr_w(1) =  (1 - MS_MeMBer_predict.r(j) + MS_MeMBer_predict.r(j) * gamma);   % first subset is allways the emtpy subset 
    dW{j}     = curr_w ;
end % end for each Bernoulli component

end % end function find_best_subsets


function [n_partitions, unique_partitions, alpha_P] = find_best_partitions(M_kk1, trellis_paths, dW, max_partitions, clutter_rate)
% Function that constructs the best partitions from the set of subsets
% returned by the find_best_subsets function.
global flag
n_partitions = 1;  % total number of partitions
partitions_list{1} = []; % we start with a single empty partition. The list of partitions will inrease after processing the subsets of each Bernoulli.
unique_partitions = zeros(1,M_kk1);  
log_w_partitions = 0; % weights of partitions, i.e. dP

for r = 1:M_kk1 % for each predicted Bernoulli component
    n_temp = 0; % local variable for the number of partitions 
    w_temp = []; 
    temp_list = {}; 
    temp_up = [];
    
    for i = 1:n_partitions % each preavious partition will branck into severall new partitions by appending subsets of the r-th Bernoulli component
        % trellis_paths{r} contains the W subsets of the r-th Bernoulli component 
        for l = 1:size(trellis_paths{r},1) 
            % check to see if an observation contained in the current
            % trellis path trellis_paths{r}(l,:) is not already contained
            % in the partition partitions_list{i}.
            if ~overlap(partitions_list{i},trellis_paths{r}(l,:))
                n_temp = n_temp + 1; % number of partitions
                temp_up(n_temp,:) = unique_partitions(i,:);
                
                if sum(trellis_paths{r}(l,:) == -1) ~=  length(trellis_paths{r}(l,:))   % if not the empty path trellis_path_index{r}(l) ~= 0                
                    % formation of a new partition by appending trellis_paths{r}(l,:)
                    % trellis_paths{r}(l,:) is a row of length s which contains the indices of observations (one for each of the s sensors) contained in the subset W 
                    temp_list{n_temp} = [partitions_list{i}; trellis_paths{r}(l,:)]; 
                    temp_up(n_temp,r) = l; % remember which subset was added to this partition  trellis_path_index{r}(l)                 
                else                    
                    % trellis_path_index{r}(l) = 0 signals that
                    % trellis_paths{r}(l,:) = [-1 -1 ... -1] (ie empty path)
                    temp_list{n_temp} = partitions_list{i};                   
                end
                w_temp(n_temp) = log_w_partitions(i) + log( dW{r}(l) ); % w_temp represents the partition weight d_P = prod(d_W)                
            end % if overlap
        end % for l each subset W associated with the Bernoulli component r
    end % for r each of the previous partitions 
    
    % keep only unique partitions
    [temp_up,index,~] = unique(temp_up,'rows');  
    temp_list = temp_list(index);
    w_temp = w_temp(index);
    n_temp = size(index,1);
    
    % keep partitions with highest weights
    n_partitions = min(max_partitions,n_temp);
    [~, index] = sort(w_temp,'descend');
    partitions_list = temp_list(index(1:n_partitions)); % update the partitions list
    unique_partitions = temp_up(index(1:n_partitions),:);
    log_w_partitions = w_temp(index(1:n_partitions));
end  % for each Bernoulli component

% compute partition weights alpha_{P}
log_alpha_P = zeros(size(log_w_partitions));
for i =1:n_partitions    
    nb_obs = sum(sum( partitions_list{i} ~= -1 )); 
    
    if flag==1
        log_alpha_P(i) = log_w_partitions(i) ;     
    else
        log_alpha_P(i) = log_w_partitions(i) + log( clutter_rate^(-nb_obs) ); 
    end    
end

 alpha_P = exp( log_alpha_P - log( sum( exp(log_alpha_P) ) ) ); % normalize


end  % end function find_best_partitions


function does_overlap = overlap(partition, path)
% function to verify if the observation subset 'path' overlaps with the 'partition'
does_overlap = 0;

p_size = size(partition,1);
nb_sensors = size(path,2);

p = 0;
while ~does_overlap && p < p_size
    p = p + 1;
    temp = partition(p,:);
    for s = 1:nb_sensors
        if (temp(s) ~= -1) && (path(s) ~= -1) && (temp(s) == path(s))
            does_overlap = 1;
            break;
        end
    end
end

end % end function overlap
 

function member_out = prune(member_in, filter)
% prune Bernoulli components with weights below a given threshold
idx = find( member_in.r > filter.comp_threshold );

member_out.M = length(idx);
member_out.r = member_in.r(idx);

member_out.m = member_in.m(:,idx);
member_out.P = member_in.P(:,:,idx);

end % end function prune


function member_out = cap(member_in, N_cap)
% Cap total number of PHD components to a specified maximum

if member_in.M > N_cap
    [~, idxsort] = sort(member_in.r,'descend');
    idxkeep = idxsort(1:N_cap);

    member_out.M = length(idxkeep);
    member_out.r = member_in.r(idxkeep);
   
    member_out.m = member_in.m(:,idxkeep);
    member_out.P = member_in.P(:,:,idxkeep);
else
    member_out = member_in;
end

end % end function cap


function [X_hat, N_hat, ind] = extract_estimates(member_in)
% State inference from the updated multi PHD distribution

% Estimate cardinality distribution
cdn_pmf = prod(1-member_in.r) * esf(member_in.r ./ (1-member_in.r));

% MAP cardinality estimate
[~, idx_max_cdn] = max(cdn_pmf);
map_cdn = idx_max_cdn-1;
N_hat = min(member_in.M, map_cdn);   % to avoid array access overflow

[~, ind] = sort(member_in.r,'descend');
ind = ind(1:N_hat);

X_hat = [];
if N_hat
    % pick highest weight components as targets
    X_hat = member_in.m(:,ind);
end

end % end function extract estimates


function display_diaginfo( k, est, filter, H_predict, H_posterior, H_prune, H_cap)
% Display diagnostic information on the fly
if ~strcmp(filter.run_flag,'silence')
    disp([' time= ',num2str(k),...
          ' #est card=' num2str(est.state_n_hat(k),4),...
          ' #comp pred=' num2str(H_predict,4),...
          ' #comp post=' num2str(H_posterior,4),...
          ' #comp pruned=',num2str(H_prune),...
          ' #comp updt=',num2str(H_cap,4)   ]);
end
end % end function display

function Xc= get_comps(X,c)

if isempty(X)
    Xc= [];
else
    Xc= X(c,:);
end

end
