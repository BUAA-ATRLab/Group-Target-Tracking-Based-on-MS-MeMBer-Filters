function [out_paths,out_w,out_m,out_P]=MDA_window(Z,s_init,s_end,MS_MeMBer_predict,model, filter, in_paths,in_w,in_m,in_P)
%初始化
out_paths=[];
out_w=[];
out_m=[];
out_P=[];

gamma = prod(1 - model.obs.pD(1:s_end));

%更改空集合的权重
curr_s=size(in_paths,2)-1;
if curr_s>0
    empty_id=repmat([MS_MeMBer_predict.M+1:2*MS_MeMBer_predict.M]',1,curr_s);
    [~,idx]=ismember(in_paths(:,2:end),empty_id,'rows');
    id=find(idx~=0);
    in_w(id)=log(MS_MeMBer_predict.r(idx(id))*prod(1-model.obs.pD(1:curr_s)));
end


%按伯努利项依次更新
% for each component
% 1:M_kk1 each bernoulli component; 
% M_kk1+1:2*M_kk1 empty measurment;
% >2*M_kk1 measurement
for i = 1:MS_MeMBer_predict.M
    idx=find(in_paths(:,1)==i);
    curr_paths=in_paths(idx,:);
    curr_w=in_w(:,idx);
    curr_m=in_m(:,idx);
    curr_P=in_P(:,:,idx);
    %传感器循环
    for s = s_init:s_end
        curr_J=length(curr_w);
        curr_Z = Z{s};
        [z_id, z_gate]= my_gate_meas(curr_Z, filter.gamma, model, MS_MeMBer_predict.m(:,i), MS_MeMBer_predict.P(:,:,i));
        z_id=z_id+2*MS_MeMBer_predict.M;
        nb_meas = size(z_gate, 2);    % current number of measurements   
        
        % Kalman updates for the for no detected observation case
        local_w = curr_w+log(1 - model.obs.pD(s));
        local_m = curr_m;
        local_P = curr_P;
        if nb_meas ~= 0   % if there is at least one measurement
          
            [qz_update, m_update, P_update] = kalman_update_multiple(z_gate, model, curr_m, curr_P);
           
%             new_w = (model.obs.pD(s)/model.obs.pdf_c/model.obs.lambda_c*repmat(curr_w, nb_meas,1) .* qz_update')';
            new_w = log(model.obs.pD(s))-log(model.obs.pdf_c)-log(model.obs.lambda_c)+repmat(curr_w, nb_meas,1)+log(qz_update');
            new_w = new_w';
            
            local_w = [local_w new_w(:)'];
            local_m = [local_m, m_update];
            local_P = cat(3, local_P, repmat(P_update, 1,1, nb_meas));  
        end
        local_paths = repmat(curr_paths, nb_meas+1, 1);
        new_obs     = kron(z_id, ones(curr_J, 1));
        new_obs     = [(i+MS_MeMBer_predict.M)*ones(curr_J,1); new_obs(:)];
        local_paths = [local_paths new_obs];      
        
        curr_w = local_w; % first association is with the void measurement
        curr_m = local_m;
        curr_P = local_P;
        curr_paths = local_paths;        
    end
    
    empty_id=repmat([MS_MeMBer_predict.M+1:2*MS_MeMBer_predict.M]',1,s_end);
    [~,idx]=ismember(curr_paths(:,2:end),empty_id,'rows');
    id=find(idx~=0);
    curr_w(id) =  log(1 - MS_MeMBer_predict.r(i) + MS_MeMBer_predict.r(i) * gamma);

    
   out_paths=[out_paths;curr_paths];
   out_w=[out_w curr_w];
   out_m = [out_m curr_m];
   out_P = cat(3, out_P, curr_P);  
end

end