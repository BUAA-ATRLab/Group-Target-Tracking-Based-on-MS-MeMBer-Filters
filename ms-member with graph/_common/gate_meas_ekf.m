function z_gate= gate_meas_ekf(z,s_indx, gamma,model,m,P)

valid_idx = [];
zlength = size(z,2);
plength = size(m,2);

for j=1:plength

        [C_ekf,U_ekf]= ekf_update_mat(model, s_indx, m(:,j));
        Sj= U_ekf*model.obs.R*U_ekf' + C_ekf*P(:,:,j)*C_ekf'; Sj= (Sj+ Sj')/2;
        
        
        Vs= chol(Sj); det_Sj= prod(diag(Vs))^2; inv_sqrt_Sj= inv(Vs);
        iSj= inv_sqrt_Sj*inv_sqrt_Sj'; 
        
        nu= z- repmat(model.obs.gen_sen_obs(s_indx,m(:,j),[]),[1 zlength]);
        
        
        dist = sum((inv_sqrt_Sj'*nu).^2);
        valid_idx = union(valid_idx,find( dist < gamma ));
        
        
        
end
z_gate = z(:,valid_idx);