function valid_idx = gate_meas_ekf_idx(z, s_idx, gamma,model,m,P)

valid_idx = [];
zlength = size(z,2);
plength = size(m,2);

eta = model.obs.gen_sen_obs(s_idx, m, []);

for j=1:plength
    [C_ekf,U_ekf]= ekf_update_mat(model,m(:,j));
    
    Sj= U_ekf*model.obs.R*U_ekf' + C_ekf*P(:,:,j)*C_ekf'; Sj= (Sj+ Sj')/2;
    Vs= chol(Sj); det_Sj= prod(diag(Vs))^2; inv_sqrt_Sj= inv(Vs);
    iSj= inv_sqrt_Sj*inv_sqrt_Sj';
    
    
    [H_ekf,U_ekf] = ekf_update_mat(model,s_indx, m);                 % user specified function for application
    S = U_ekf*model.obs.R*U_ekf' + H_ekf*P*H_ekf'; S = (S+ S')/2;
    Vs= chol(S); det_S = prod(diag(Vs))^2; inv_sqrt_S = inv(Vs); iS = inv_sqrt_S*inv_sqrt_S';
    
    K  = P*H_ekf'*iS;
    
    innov = model.obs.angle_wrap( z-repmat(eta,[1 size(z,2)]) );
    
    
    dist= sum((inv_sqrt_Sj'*innov).^2);
    valid_idx= unique_faster([ valid_idx find( dist < gamma )]);
end
valid_idx = valid_idx(:)';

end