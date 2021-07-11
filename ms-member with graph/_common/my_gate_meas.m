function [z_id, z_gate]= my_gate_meas(z, gamma, model, m, P)

z_id = [];
z_gate=[];
zlength = size(z,2);



Sj = model.obs.R1 + model.obs.H*P*model.obs.H';
Vs = chol(Sj); det_Sj= prod(diag(Vs))^2;  inv_sqrt_Sj= inv(Vs);
iSj = inv_sqrt_Sj*inv_sqrt_Sj';
nu = z- model.obs.H*repmat(m,[1 zlength]);
dist = sum((inv_sqrt_Sj'*nu).^2);
z_id = union(z_id,find( dist < gamma ));
z_gate = z(:,z_id);


end