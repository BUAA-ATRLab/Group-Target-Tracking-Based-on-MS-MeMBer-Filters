function z_gate= gate_meas_ukf(z, s_indx, gamma, model, m, P, alpha, kappa, beta)

valid_idx = [];
zlength = size(z,2);
plength = size(m,2);

for j=1:plength
        [X_ukf,u] = ut( [m(:,j); zeros(model.obs.z_dim,1) ], blkdiag(P(:,:,j),model.obs.R), alpha, kappa );
        Z_pred = model.obs.gen_sen_obs( s_indx, X_ukf(1:model.kin.x_dim,:), X_ukf((model.kin.x_dim+1):(model.kin.x_dim+model.obs.z_dim),:) );
        eta = Z_pred*u(:); Sj_temp = Z_pred- repmat(eta,[1 length(u)]); u(1)= u(1)+(1-alpha^2+beta);
        Sj = Sj_temp*diag(u)*Sj_temp';
        
        Vs = chol(Sj); det_Sj = prod(diag(Vs))^2; inv_sqrt_Sj= inv(Vs);
        iSj = inv_sqrt_Sj*inv_sqrt_Sj'; 
        nu = z - repmat(model.obs.gen_sen_obs( s_indx, m(:,j), []),[1 zlength]);
        dist = sum((inv_sqrt_Sj'*nu).^2);
        
        valid_idx = union(valid_idx,find( dist < gamma ));
end
z_gate = z(:,valid_idx);

end