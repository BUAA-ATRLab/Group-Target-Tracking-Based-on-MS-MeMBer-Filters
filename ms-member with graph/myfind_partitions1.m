% Grasp最大权重独立集算法求最优划分P
function out = myfind_partitions1(model, filter, MS_MeMBer_predict, Z)    
%% Part 1 find likely assignment for each Bernoulli component
M_kk1 = MS_MeMBer_predict.M; % number of predicted Bernoulli components
gamma = prod(1 - model.obs.pD);

% initialize output variables
best_paths = cell(1,M_kk1);
m_paths = cell(1,M_kk1);
P_paths = cell(1,M_kk1);
dW = cell(1,M_kk1);

% for each component
% 1:M_kk1 each bernoulli component; 
% M_kk1+1:2*M_kk1 empty measurment;
% >2*M_kk1 measurement
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
        [z_id, z_gate]= my_gate_meas(curr_Z, filter.gamma, model, MS_MeMBer_predict.m(:,j), MS_MeMBer_predict.P(:,:,j));
        z_id=z_id+2*M_kk1;
        nb_meas = size(z_gate, 2);    % current number of measurements   
        
        % Kalman updates for the for no detected observation case
        local_w = curr_w*(1 - model.obs.pD(s));
        local_m = curr_m;
        local_P = curr_P;
        
        if nb_meas ~= 0   % if there is at least one measurement
          
            [qz_update, m_update, P_update] = kalman_update_multiple(z_gate, model, curr_m, curr_P);
           
            new_w = (model.obs.pD(s)/model.obs.pdf_c/model.obs.lambda_c*repmat(curr_w, nb_meas,1) .* qz_update')';
            
            local_w = [local_w new_w(:)'];
            local_m = [local_m, m_update];
            local_P = cat(3, local_P, repmat(P_update, 1,1, nb_meas));
        end
 
        local_paths = repmat(curr_paths, nb_meas+1, 1);
        new_obs     = kron(z_id, ones(curr_J, 1));
        new_obs     = [(j+M_kk1)*ones(curr_J,1); new_obs(:)];
        local_paths = [local_paths new_obs];
        

        % select best observation sequences till now 
        % always include the empty observation path first
        [~,index] = sort(local_w(2:end),'descend');
        % select only the best W_max measurements as candidates for the j component
%         J = min(length(index), filter.W_max);  % update number of components propagated
        J = length(index);
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
    dW{j}     = log(curr_w) ;
end % end for each Bernoulli component
%% Part 2 S-D assignment
sub=[];
valCost=[];
m_all=[];
P_all=[];
for j = 1:M_kk1
    temp=[j*ones(size(best_paths{j},1),1) best_paths{j}];
    sub=[sub;temp];
    valCost=[valCost;dW{j}(:)];
    m_all = [m_all m_paths{j}];
    P_all = cat(3, P_all, P_paths{j});
end
minvalue=min(min(valCost));
if isinf(minvalue)
    minvalue=-3000;
end
val=valCost-minvalue+1;
[gam1, cost1] = k_best_myGrasp(sub, val, filter.SD_max, filter.sGrasp);
% [gam, cost] = k_best_assignment(sub, -val, filter.SD_max);
log_alpha_P=zeros(1,size(gam1,3));
for j=1:size(gam1,3);
    [~,Locb] = ismember(gam1(:,:,j),sub,'rows');
    % transform assignment mat
    gam_temp=gam1(:,2:end,j);
    idx0=find(gam_temp>M_kk1&gam_temp<=2*M_kk1);
    gam_temp(idx0)=0;
    idx1=find(gam_temp>2*M_kk1);
    gam_temp(idx1)=gam_temp(idx1)-2*M_kk1;
    gam1(:,2:end,j)=gam_temp;
    %% out
    out.sub{j}=gam1(:,:,j);
    out.valCost{j}=valCost(Locb,:);
    out.m{j}=m_all(:,Locb);
    out.P{j}=P_all(:,:,Locb);
    log_alpha_P(j)=sum(valCost(Locb,:));
end
out.alpha_P = exp( log_alpha_P - log( sum( exp(log_alpha_P) ) ) ); % normalize
end